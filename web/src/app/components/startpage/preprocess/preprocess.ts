import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { Component, inject, signal } from '@angular/core';
import { ApiService, ProgressEvent} from '../../../shared/api.service';
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
  isProcessing = signal<boolean>(false);
  processedFileName = signal('');

  currentJobId: string | null = null;
  progress = 0;
  currentStage = '';
  currentMessage = '';
  logs: ProgressEvent[] = [];

  constant5Region: string = '';
  constant3Region: string = '';
  sequenceLengthMin: number = 10;
  sequenceLengthMax: number = 100;
  sequenceLengthRange: number = 500;
  maxAllowedError: number = 0.005;
  maxErrorRange: number = 1;
  progressSubscription: any;

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
    this.logs = [];
    this.progress = 0;

    const params = {
      input_path: this.savedFileName,
      const5p: this.constant5Region,
      const3p: this.constant3Region,
      min_length: this.sequenceLengthMin,
      max_length: this.sequenceLengthMax,
      max_error: this.maxAllowedError,
      output_format: this.downloadFormat
    };

    this.apiService.startPreprocessing(params).subscribe({
      next: (jobId) => {
        this.currentJobId = jobId;
        this.subscribeToProgress(jobId);
      },
      error: (error) => {
        console.error('Failed to start preprocessing:', error);
        this.isProcessing.set(false);
        alert('Failed to start preprocessing');
      }
    });
  }

  private subscribeToProgress(jobId: string): void {
    this.progressSubscription = this.apiService
      .subscribeToProgress(jobId)
      .subscribe({
        next: (event: ProgressEvent) => {
          if (event.stage === 'complete') {
            this.isProcessing.set(false);
            this.processedFileName.set(event.data?.output_path ?? '');
          } else if (event.stage === 'error') {
            this.isProcessing.set(false);
            console.error('SSE error:', event.message);
          } else{
            this.currentStage = event.stage ?? '';
            this.currentMessage = event.message ?? '';
            if (event.progress != null) this.progress = event.progress; // 0..1
          } 
          this.logs.push(event);
          setTimeout(() => this.scrollToBottom(), 50);
        },
        error: (error) => {
          console.error('Progress stream error:', error);
          this.isProcessing.set(false);
          alert('Connection to server lost');
        },
        complete: () => {
          this.isProcessing.set(false);
          console.log('Processing completed');
        }
      });
  }

  cancelProcessing(): void {
    if (this.currentJobId) {
      this.apiService.cancelJob(this.currentJobId);
      this.progressSubscription?.unsubscribe();
      this.isProcessing.set(false);
      this.currentJobId = null;
    }
  }

  private scrollToBottom(): void {
    const logContainer = document.getElementById('log-container');
    if (logContainer) {
      logContainer.scrollTop = logContainer.scrollHeight;
    }
  }

  getStageIcon(stage: string): string {
    const icons: { [key: string]: string } = {
      'init': '🚀',
      'cutadapt': '✂️',
      'qc': '🔍',
      'write': '💾',
      'convert': '🔄',
      'complete': '✅',
      'error': '❌'
    };
    return icons[stage] || '⏳';
  }

  ngOnDestroy(): void {
    this.cancelProcessing();
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
