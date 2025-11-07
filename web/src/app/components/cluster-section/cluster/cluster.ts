import { Component, inject, signal } from '@angular/core';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';

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
      
      this.apiService.clusterLed(params).subscribe({
        next: (response) => {
          this.isProcessing.set(false);
          this.processedFileName.set(response.result);
          console.log('Clustering completed:', response.result);
        },
        error: (error) => {
          console.error('Clustering failed:', error);
          this.isProcessing.set(false);
          
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
        }
      });
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
