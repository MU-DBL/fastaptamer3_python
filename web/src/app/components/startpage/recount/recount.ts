import { Component, inject, signal, output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';

@Component({
  selector: 'app-recount',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './recount.html',
  styleUrl: './recount.scss'
})
export class Recount {
  private apiService = inject(ApiService);

  // Output events to emit to parent component
  resultsReady = output<any[]>();
  showReadsPerRankModal = output<{ data: any[], params: any }>();
  showSeqLengthModal = output<{ data: any, params: any }>();

  // File upload states for two FASTA files
  selectedFile1: File | null = null;
  savedFileName1: string = '';
  uploadComplete1: boolean = false;

  selectedFile2: File | null = null;
  savedFileName2: string = '';
  uploadComplete2: boolean = false;

  // Form values
  normalizeValue: string = '1e+06';
  downloadFormat: string = 'fasta';

  // Use signals for reactive state
  isProcessing = signal(false);
  processedFileName = signal('');

  // Table data
  private tableData: any[] = [];

  // Plot parameters
  minReadsToPlot: number = 10;
  maxRankToPlot: number = 100;

  // Reads per Rank plot customization
  adjustReadsPerRank: string = 'no';
  rprXAxis: string = 'Ranks of unique sequences';
  rprYAxis: string = 'Total reads per unique sequence';
  rprTitle: string = 'Read count for each rank';
  rprLineColor: string = '#87CEEB';

  // Sequence-length histogram customization
  adjustSeqLengthHistogram: string = 'no';
  histXAxis: string = 'Sequence length';
  histYAxis1: string = 'Unique sequences';
  histYAxis2: string = 'Read count';
  histTitle: string = 'Sequence-length histogram';
  histBarOutline: string = '#000000';
  histBarFill: string = '#87CEEB';
  histBarFill2: string = '#FFA500';

  // Statistics
  totalSequences = signal(0);
  uniqueSequences = signal(0);

  // Helper method for slider label formatting
  formatLabel(value: number): string {
    return `${value}`;
  }

  // File 1 upload handlers
  onFile1Selected(result: FileUploadResult): void {
    this.selectedFile1 = result.file;
    this.processedFileName.set('');
    this.tableData = [];
    console.log('File 1 selected:', result.fileName);
  }

  onUpload1Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete1 = true;
      this.savedFileName1 = result.savedFileName;
      console.log('Upload 1 complete:', result.savedFileName);
    } else if (result.error) {
      console.error('Upload 1 failed:', result.error);
    }
  }

  // File 2 upload handlers
  onFile2Selected(result: FileUploadResult): void {
    this.selectedFile2 = result.file;
    this.processedFileName.set('');
    this.tableData = [];
    console.log('File 2 selected:', result.fileName);
  }

  onUpload2Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete2 = true;
      this.savedFileName2 = result.savedFileName;
      console.log('Upload 2 complete:', result.savedFileName);
    } else if (result.error) {
      console.error('Upload 2 failed:', result.error);
    }
  }

  onStart(): void {
    if (!this.uploadComplete1 || !this.savedFileName1) {
      console.warn('Please upload the first FASTA file!');
      return;
    }

    if (!this.uploadComplete2 || !this.savedFileName2) {
      console.warn('Please upload the second FASTA file!');
      return;
    }

    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.tableData = [];

    const params = {
      input_path_1: this.savedFileName1,
      input_path_2: this.savedFileName2,
      scaling_factor: parseFloat(this.normalizeValue),
      output_format: this.downloadFormat
    };

    console.log('Starting recount with parameters:', params);

    // Chain the operations using RxJS operators
    this.apiService.recount(params).pipe(
      tap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Recount complete:', response.result);
        }
      }),
      switchMap(response => {
        // Automatically load results after successful recount
        if (response.status === 'ok' && response.result) {
          return this.apiService.downloadFile(response.result).pipe(
            tap(blob => this.parseFileBlob(blob, response.result))
          );
        }
        return of(null);
      }),
      catchError(error => {
        const errorMsg = error.error?.detail || 'Recount failed';
        console.error('Recount error:', errorMsg);
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
      const headers = lines[0].split(',');

      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',');
        if (values.length >= 6) {
          this.tableData.push({
            id: values[0],
            rank: parseInt(values[1]),
            reads: parseInt(values[2]),
            rpm: parseInt(values[3]),
            length: parseInt(values[4]),
            seqs: values[5]
          });
        }
      }
    } else if (isFasta) {
      // Parse FASTA format
      const lines = content.split('\n').filter(line => line.trim());
      let currentId = '';
      let currentSeq = '';

      for (const line of lines) {
        if (line.startsWith('>')) {
          if (currentId && currentSeq) {
            this.addFastaEntry(currentId, currentSeq);
          }
          currentId = line.substring(1);
          currentSeq = '';
        } else {
          currentSeq += line.trim();
        }
      }

      // Add the last entry
      if (currentId && currentSeq) {
        this.addFastaEntry(currentId, currentSeq);
      }
    }

    // Update statistics
    this.uniqueSequences.set(this.tableData.length);
    const totalReads = this.tableData.reduce((sum, item) => sum + item.reads, 0);
    this.totalSequences.set(totalReads);

    // Emit the parsed data to the parent component
    this.resultsReady.emit(this.tableData);
  }

  addFastaEntry(id: string, sequence: string): void {
    // Parse ID format: "rank=1;read=20;RPU=200000"
    const parts = id.split(';');
    const rank = parts[0]?.split('=')[1] || '';
    const reads = parts[1]?.split('=')[1] || '';
    const rpm = parts[2]?.split('=')[1] || '';

    this.tableData.push({
      id: id,
      rank: parseInt(rank),
      reads: parseInt(reads),
      rpm: parseInt(rpm),
      length: sequence.length,
      seqs: sequence
    });
  }

  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) {
      console.warn('No file available for download. Please run recount first.');
      return;
    }

    console.log('Downloading file:', filename);

    this.apiService.downloadFile(filename).subscribe({
      next: (blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        window.URL.revokeObjectURL(url);
        console.log('Download started');
      },
      error: (error) => {
        const errorMsg = error.error?.detail || 'Download failed';
        console.error('Download error:', errorMsg);
      }
    });
  }

  // Plot methods - emit events to parent component
  openReadsPerRankPlot(): void {
    if (this.tableData.length === 0) {
      alert('No data available for plotting. Please run the recount process first.');
      return;
    }

    // Validate slider values
    if (this.minReadsToPlot < 0) {
      alert('Minimum number of reads cannot be negative.');
      return;
    }

    if (this.maxRankToPlot < this.minReadsToPlot) {
      alert('Maximum rank must be greater than or equal to minimum reads.');
      return;
    }

    if (this.maxRankToPlot < 1) {
      alert('Maximum rank must be at least 1.');
      return;
    }

    // Check if max rank exceeds available data
    const maxAvailableRank = Math.max(...this.tableData.map(item => item.rank));
    const data = this.getReadsPerRankData();

    if (data.length === 0) {
      alert('No data points match the specified criteria. Please adjust the min reads or max rank values.');
      return;
    }

    // Inform user if max rank was adjusted
    if (this.maxRankToPlot > maxAvailableRank) {
      console.info(`Max rank adjusted from ${this.maxRankToPlot} to ${maxAvailableRank} (maximum available rank)`);
    }

    const params = {
      title: this.rprTitle,
      xAxisLabel: this.rprXAxis,
      yAxisLabel: this.rprYAxis,
      lineColor: this.rprLineColor
    };

    this.showReadsPerRankModal.emit({ data, params });
  }

  openSeqLengthHistogram(): void {
    if (this.tableData.length === 0) {
      alert('No data available for plotting. Please run the recount process first.');
      return;
    }

    const data = this.getSequenceLengthData();

    if (data.unique.length === 0 && data.total.length === 0) {
      alert('No sequence length data available for plotting.');
      return;
    }

    const params = {
      title: this.histTitle,
      xAxisLabel: this.histXAxis,
      yAxis1Label: this.histYAxis1,
      yAxis2Label: this.histYAxis2,
      barOutline: this.histBarOutline,
      barFill: this.histBarFill,
      barFill2: this.histBarFill2
    };

    this.showSeqLengthModal.emit({ data, params });
  }

  getReadsPerRankData(): any[] {
    // Get all data within the rank range, regardless of minReadsToPlot
    const maxAvailableRank = Math.max(...this.tableData.map(item => item.rank));
    const effectiveMaxRank = Math.min(this.maxRankToPlot, maxAvailableRank);

    return this.tableData
      .filter(item => item.reads >= this.minReadsToPlot && item.rank <= effectiveMaxRank)
      .map(item => ({ rank: item.rank, reads: item.reads }))
      .sort((a, b) => a.rank - b.rank);
  }

  getSequenceLengthData(): any {
    const uniqueLengths: { [key: number]: number } = {};
    const totalReadsByLength: { [key: number]: number } = {};

    this.tableData.forEach(item => {
      uniqueLengths[item.length] = (uniqueLengths[item.length] || 0) + 1;
      totalReadsByLength[item.length] = (totalReadsByLength[item.length] || 0) + item.reads;
    });

    return {
      unique: Object.keys(uniqueLengths).map(len => ({ length: parseInt(len), count: uniqueLengths[parseInt(len)] })),
      total: Object.keys(totalReadsByLength).map(len => ({ length: parseInt(len), count: totalReadsByLength[parseInt(len)] }))
    };
  }
}
