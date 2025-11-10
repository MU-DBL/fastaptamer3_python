import { Component, inject, signal, output } from '@angular/core';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';

export interface ClusterResultsEvent {
  data: any[];
  inputFile?: string;
}

@Component({
  selector: 'app-cluster',
  imports: [   
    CommonModule,
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster.html',
  styleUrl: './cluster.scss',
})
export class Cluster {
    private apiService = inject(ApiService);

    // Output event for results ready
    resultsReady = output<ClusterResultsEvent>();

    downloadFormat: string = 'fasta';
    keepNonClusteredSequence: string = 'no';

    isProcessing = signal(false);
    processedFileName = signal('');
    
    selectedFile: File | null = null;
    savedFileName: string = '';
    uploadComplete: boolean = false;

    number_reads_to_cluster: number = 10;
    led: number = 7;
    number_of_cluster: number = 20;

    onFileSelected(result: FileUploadResult): void {
      this.selectedFile = result.file;
      this.processedFileName.set('');
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

      this.isProcessing.set(true);
      this.processedFileName.set('');

      const params = {
        input_path: this.savedFileName,
        output_format: this.downloadFormat,
        min_reads: this.number_reads_to_cluster,
        max_led: this.led,
        total_clusters: this.number_of_cluster,
        keep_nc: this.keepNonClusteredSequence === 'yes'
      };

      console.log('Starting clustering with params:', params);
      
      this.apiService.clusterLed(params).pipe(
        tap(response => {
          if (response.status === 'ok' && response.result) {
            this.processedFileName.set(response.result);
            console.log('Clustering completed:', response.result);
          }
        }),
        switchMap(response => {
          // Automatically load and parse results after successful clustering
          if (response.status === 'ok' && response.result) {
            return this.apiService.downloadFile(response.result).pipe(
              tap(blob => this.parseClusterFile(blob, response.result))
            );
          }
          return of(null);
        }),
        catchError(error => {
          console.error('Clustering failed:', error);
          
          let errorMessage = 'Clustering failed: ';
          if (error.error?.detail) {
            errorMessage += error.error.detail;
          } else {
            errorMessage += 'Unknown error';
          }
          
          // Show user-friendly message for common errors
          if (error.error?.detail?.includes('No sequences remain after filtering')) {
            errorMessage += '\n\nTip: Try reducing the "Min. number of reads to cluster" value.';
          }
          
          alert(errorMessage);
          return of(null);
        }),
        finalize(() => {
          this.isProcessing.set(false);
        })
      ).subscribe();
    }

    parseClusterFile(blob: Blob, filename: string): void {
      const reader = new FileReader();
      reader.onload = (e: any) => {
        const text = e.target.result;
        this.parseResultFile(text, filename);
      };
      reader.readAsText(blob);
    }

    parseResultFile(content: string, filename: string): void {
      const isCsv = filename.endsWith('.csv');
      const clusterData: any[] = [];
      
      if (isCsv) {
        // Parse CSV content (similar to diversity component but for cluster results)
        const lines = content.split('\n').filter(line => line.trim());
        if (lines.length > 1) {
          const headers = lines[0].split(',');
          
          for (let i = 1; i < lines.length; i++) {
            const values = lines[i].split(',');
            if (values.length >= 4) {
              const rowData = {
                sequenceID: values[0] || '',
                cluster: parseInt(values[1]) || 0,
                readCount: parseInt(values[2]) || 0,
                sequence: values[3] || '',
                // Add additional fields if they exist
                ...(values.length > 4 && { additionalData: values.slice(4) })
              };
              clusterData.push(rowData);
            }
          }
        }
      } else if (filename.endsWith('.fasta') || filename.endsWith('.fa')) {
        // Parse FASTA content
        const fastaEntries = this.parseFasta(content);
        fastaEntries.forEach((entry, index) => {
          // Extract cluster information from header if available
          const clusterMatch = entry.header.match(/Cluster[_\s]*(\d+)/i);
          const readCountMatch = entry.header.match(/reads[_\s]*(\d+)/i);
          
          clusterData.push({
            sequenceID: entry.header,
            cluster: clusterMatch ? parseInt(clusterMatch[1]) : index + 1,
            readCount: readCountMatch ? parseInt(readCountMatch[1]) : 1,
            sequence: entry.sequence
          });
        });
      }

      // Emit results to parent component for display in right section
      this.resultsReady.emit({
        data: clusterData,
        inputFile: this.savedFileName
      });
      
      console.log('Parsed cluster data:', clusterData.length, 'sequences');
    }

    private parseFasta(content: string): { header: string; sequence: string }[] {
      const entries: { header: string; sequence: string }[] = [];
      const lines = content.split('\n');
      let currentEntry: { header: string; sequence: string } | null = null;
      
      for (const line of lines) {
        const trimmedLine = line.trim();
        if (trimmedLine.startsWith('>')) {
          // New sequence header
          if (currentEntry) {
            entries.push(currentEntry);
          }
          currentEntry = {
            header: trimmedLine.substring(1),
            sequence: ''
          };
        } else if (currentEntry && trimmedLine) {
          // Add to current sequence
          currentEntry.sequence += trimmedLine;
        }
      }
      
      // Don't forget the last entry
      if (currentEntry) {
        entries.push(currentEntry);
      }
      
      return entries;
    }

    onDownload(): void {
      const filename = this.processedFileName();
      if (!filename) {
        console.warn('No file available for download. Please run clustering first.');
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
          alert('Download failed: ' + errorMsg);
        }
      });
    }
}
