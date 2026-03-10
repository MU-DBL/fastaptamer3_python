import { Component, inject, signal, ChangeDetectorRef, OnDestroy } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';
import { Table, TableConfig } from '../../common/table/table';

@Component({
  selector: 'app-motif-omit',
  imports: [ 
    CommonModule,
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './motif-omit.html',
  styleUrl: './motif-omit.scss'
})
export class MotifOmit implements OnDestroy {

  tableConfig: TableConfig = {
    columns: [
      { key: 'id', label: 'id' },
      { key: 'rank', label: 'Rank' },
      { key: 'reads', label: 'Reads' },
      { key: 'rpm', label: 'RPM' },
      { key: 'seqs', label: 'Sequence' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  private apiService = inject(ApiService);
  private cdr = inject(ChangeDetectorRef);

  selectedFile: File | null = null;
  savedFileName: string = '';
  fileName: string = 'FASTA file';
  motifPattern: string = '';
  partialMatch: string = 'no';
  motifType: string = 'Nucleotide';
  downloadFormat: string = 'fasta';
  uploadComplete: boolean = false;
  
  // Use signals for reactive state that affects the template
  isProcessing = signal(false);
  processedFileName = signal('');

  // Table data
  tableData: any[] = [];

  ngOnDestroy(): void {
    if (this.savedFileName) {
      this.apiService.deleteFile(this.savedFileName).subscribe();
    }
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
  }

  onFileSelected(result: FileUploadResult): void {
    if (this.savedFileName) {
      this.apiService.deleteFile(this.savedFileName).subscribe();
    }
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
    this.selectedFile = result.file;
    this.savedFileName = '';
    this.processedFileName.set('');
    this.tableData = [];
    console.log('File selected:', result.fileName);
  }

  onUploadComplete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete = true;
      this.savedFileName = result.savedFileName;
      console.log('Upload complete:', result.savedFileName);
    } else if (result.error) {
      console.error('Upload failed:', result.error);
    }
  }

  onStart(): void {
    if (!this.uploadComplete || !this.savedFileName) {
      console.warn('Please upload a file first!');
      return;
    }

    if (!this.motifPattern || this.motifPattern.trim() === '') {
      console.warn('Please enter a motif pattern!');
      alert('Please enter at least one motif pattern!');
      return;
    }

    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.tableData = [];

    const params = {
      input_path: this.savedFileName,
      motif: this.motifPattern.trim(),
      partial: this.partialMatch === 'yes',
      motif_type: this.motifType,
      output_format: this.downloadFormat
    };

    console.log('Starting motif omit with parameters:', params);

    // Chain the operations using RxJS operators
    this.apiService.motifOmit(params).pipe(
      tap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Motif omit complete:', response.result);
        }
      }),
      switchMap(response => {
        // Automatically load results after successful omit
        if (response.status === 'ok' && response.result) {
          return this.apiService.downloadFile(response.result).pipe(
            tap(blob => this.parseFileBlob(blob, response.result))
          );
        }
        return of(null);
      }),
      catchError(error => {
        const errorMsg = error.error?.detail || 'Motif omit failed';
        console.error('Motif omit error:', errorMsg);
        alert(`Error: ${errorMsg}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  parseFileBlob(blob: Blob, filename: string): void {
    const reader = new FileReader();
    reader.onload = (e: any) => {
      const text = e.target.result;
      this.parseResultFile(text, filename);
      // Manually trigger change detection since FileReader runs outside Angular's zone
      this.cdr.detectChanges();
    };
    reader.readAsText(blob);
  }

  parseResultFile(content: string, filename: string): void {
    const isFasta = filename.endsWith('.fasta') || filename.endsWith('.fa');
    const isCsv = filename.endsWith('.csv');
    
    this.tableData = [];
    
    if (isCsv) {
      // Parse CSV
      const lines = content.split('\n').filter(line => line.trim());
      const headers = lines[0].split(',').map(h => h.trim());
      
      // Find column indices
      const idIdx = headers.findIndex(h => h === 'ID');
      const rankIdx = headers.findIndex(h => h === 'Rank');
      const readsIdx = headers.findIndex(h => h === 'Reads');
      const rpuIdx = headers.findIndex(h => h === 'RPU');
      const seqIdx = headers.findIndex(h => h === 'sequences');
      
      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',');
        if (values.length > 0 && idIdx >= 0) {
          this.tableData.push({
            id: values[idIdx] || '',
            rank: rankIdx >= 0 ? parseInt(values[rankIdx]) : 0,
            reads: readsIdx >= 0 ? parseInt(values[readsIdx]) : 0,
            rpm: rpuIdx >= 0 ? parseInt(values[rpuIdx]) : 0,
            seqs: seqIdx >= 0 ? values[seqIdx] : ''
          });
        }
      }
    } else if (isFasta) {
      // Parse FASTA
      const lines = content.split('\n').filter(line => line.trim() !== '');
      for (let i = 0; i < lines.length; i += 2) {
        if (lines[i] && lines[i].startsWith('>') && lines[i + 1]) {
          const header = lines[i].substring(1).trim();
          const sequence = lines[i + 1].trim();
          
          // Parse header (format: Rank=1;Reads=100;RPU=1000)
          const parts: any = {};
          header.split(';').forEach(part => {
            const [key, value] = part.split('=');
            if (key && value) {
              parts[key.trim()] = value.trim();
            }
          });
          
          this.tableData.push({
            id: header,
            rank: parseInt(parts['Rank'] || '0'),
            reads: parseInt(parts['Reads'] || '0'),
            rpm: parseInt(parts['RPU'] || '0'),
            seqs: sequence
          });
        }
      }
    }
    
    console.log('Parsed results:', this.tableData.length, 'sequences');
  }

  onDownload(): void {
    if (!this.processedFileName()) {
      console.warn('No processed file available for download');
      return;
    }

    this.apiService.downloadFile(this.processedFileName()).subscribe({
      next: (blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = this.processedFileName();
        link.click();
        window.URL.revokeObjectURL(url);
        console.log('File downloaded:', this.processedFileName());
      },
      error: (error) => {
        console.error('Download error:', error);
        alert('Failed to download file');
      }
    });
  }
}
