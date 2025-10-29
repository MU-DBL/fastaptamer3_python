import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { Component, inject, signal } from '@angular/core';
import { ApiService } from '../../../shared/api.service';
import { finalize } from 'rxjs/operators';

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
  
  // Use signals for reactive state
  isProcessing = signal(false);
  processedFileName = signal('');

  constant5Region: string = '';
  constant3Region: string = '';
  sequenceLengthMin: number = 10;
  sequenceLengthMax: number = 100;
  sequenceLengthRange: number = 500;
  maxAllowedError: number = 0.005;
  maxErrorRange: number = 1;

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
      const5p: this.constant5Region,
      const3p: this.constant3Region,
      min_length: this.sequenceLengthMin,
      max_length: this.sequenceLengthMax,
      max_error: this.maxAllowedError,
      output_format: this.downloadFormat
    };

    console.log('Starting preprocessing with parameters:', params);

    this.apiService.preprocess(params).pipe(
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe({
      next: (response) => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Preprocessing complete:', response.result);
        }
      },
      error: (error) => {
        const errorMsg = error.error?.detail || 'Preprocessing failed';
        console.error('Preprocessing error:', errorMsg);
      }
    });
  }

  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) {
      console.warn('No file available for download. Please run preprocessing first.');
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

  formatLabel(value: number): string {
    return `${value}`;
  }

  formatErrorLabel(value: number): string {
    return value.toFixed(3);
  }
}
