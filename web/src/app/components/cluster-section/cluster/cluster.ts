import { Component, inject, signal, output, ChangeDetectorRef } from '@angular/core';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';
import { Table ,TableConfig} from '../../common/table/table';

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
    Table,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster.html',
  styleUrl: './cluster.scss',
})

export class Cluster {

    private cdr = inject(ChangeDetectorRef);
    private apiService = inject(ApiService);

    tableConfig: TableConfig = {
        columns: [
          { key: 'sequenceID', label: 'ID', width:"150px"},               
          { key: 'cluster', label: 'Cluster', exect_match:true },  
          { key: 'rankInCluster', label: 'Rank In Cluster', exect_match:true}, 
          { key: 'led', label: 'LED' },          
          { key: 'readCount', label: 'Reads' },             
          { key: 'rank', label: 'Rank' },                   
          { key: 'rpu', label: 'RPU' },                     
          { key: 'sequence', label: 'Sequence' }                         
        ],
        initialPageSize: 10,
        pageSizeOptions: [10, 25, 50, 100]
      };

    clusterData: any[] = [];

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
      this.clusterData = []
      const reader = new FileReader();
      reader.onload = (e: any) => {
        const text = e.target.result;
        this.parseResultFile(text, filename);
      };
      reader.readAsText(blob);
    }

  parseResultFile(content: string, filename: string): void {
  // 1. Create a temporary local array
    const newData: any[] = []; 

    const isCsv = filename.endsWith('.csv');

    if (isCsv) {
      const lines = content.split('\n').filter(line => line.trim());
      if (lines.length > 1) {
        // Skip header (i=1)
        for (let i = 1; i < lines.length; i++) {
          const values = lines[i].split(',');
          if (values.length >= 4) {
            const rowData = {
              sequenceID: values[0] || '',
              cluster: parseInt(values[1]) || 0,
              readCount: parseInt(values[2]) || 0,
              sequence: values[3] || '',
              ...(values.length > 4 && { additionalData: values.slice(4) })
            };
            // 2. Push to temporary array
            newData.push(rowData); 
          }
        }
      }
    } else if (filename.endsWith('.fasta') || filename.endsWith('.fa')) {
      const fastaEntries = this.parseFasta(content);
      
      fastaEntries.forEach((entry) => {
        const getParam = (key: string, isFloat = false) => {
          const regex = new RegExp(`${key}[=_\\s]*([0-9.]+)`, 'i');
          const match = entry.header.match(regex);
          return (match && match[1]) ? (isFloat ? parseFloat(match[1]) : parseInt(match[1])) : 0;
        };

        // 2. Push to temporary array
        newData.push({
          sequenceID: entry.header,
          sequence: entry.sequence,
          rank: getParam('Rank'),
          readCount: getParam('Reads'),
          rpu: getParam('RPU', true),
          cluster: getParam('Cluster'),
          rankInCluster: getParam('RankInCluster'),
          led: getParam('LED')
        });
      });
      }
      // console.log(newData)
    this.clusterData = newData; 
    this.cdr.detectChanges();
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
