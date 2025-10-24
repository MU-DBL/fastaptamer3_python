import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { Component, inject } from '@angular/core';
import { ApiService } from '../../../shared/api.service';

@Component({
  selector: 'app-preprocess',
  imports: [  
    CommonModule, 
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS],
  templateUrl: './preprocess.html',
  styleUrl: './preprocess.scss'
})

export class Preprocess {
  private apiService = inject(ApiService);

  selectedFile: File | null = null;
  savedFileName: string = '';
  fileName: string = 'FASTQ or FASTA file';
  downloadFormat: string = 'fasta';
  uploadComplete: boolean = false;
  isUploading: boolean = false;
  isProcessing: boolean = false;
  processedFileName: string = '';

  constant5Region: string = '';
  constant3Region: string = '';
  sequenceLengthMin: number = 10;
  sequenceLengthMax: number = 100;
  sequenceLengthRange: number = 500;
  maxAllowedError: number = 0.005;
  maxErrorRange: number = 1;

  onFileSelected(result: FileUploadResult): void {
    this.selectedFile = result.file;
    this.processedFileName = '';
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

    this.isProcessing = true;
    this.processedFileName = '';

    const params = {
      input_path: this.savedFileName,
      const5p: this.constant5Region,
      const3p: this.constant3Region,
      min_length: this.sequenceLengthMin,
      max_length: this.sequenceLengthMax,
      max_error: this.maxAllowedError,
      output_format: this.downloadFormat
    };

    console.log('Starting preprocessing with parameters:', params);

    this.apiService.preprocess(params).subscribe({
      next: (response) => {
        this.isProcessing = false;
        if (response.status === 'ok' && response.result) {
          this.processedFileName = response.result;
          console.log('Preprocessing complete:', response.result);
        }
      },
      error: (error) => {
        this.isProcessing = false;
        const errorMsg = error.error?.detail || 'Preprocessing failed';
        console.error('Preprocessing error:', errorMsg);
      }
    });
  }

  onDownload(): void {
    if (!this.processedFileName) {
      console.warn('No file available for download. Please run preprocessing first.');
      return;
    }

    console.log('Downloading file:', this.processedFileName);

    this.apiService.downloadFile(this.processedFileName).subscribe({
      next: (blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = this.processedFileName;
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

  formatLabel(value: number): string {
    return `${value}`;
  }

  formatErrorLabel(value: number): string {
    return value.toFixed(3);
  }
}
