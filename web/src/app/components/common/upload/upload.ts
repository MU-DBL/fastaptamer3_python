import { Component, EventEmitter, Input, Output, inject } from '@angular/core';
import { CommonModule } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';

export interface FileUploadResult {
  file: File;
  fileName: string;
  savedFileName?: string;
  uploadComplete: boolean;
  error?: string;
}

@Component({
  selector: 'app-upload',
  imports: [
    CommonModule,
    ...MATERIAL_IMPORTS],
  templateUrl: './upload.html',
  styleUrl: './upload.scss'
})
export class Upload {
  @Input() acceptedFileTypes: string = '.fasta,.fastq,.fa,.fq';
  @Input() buttonText: string = 'Browse...';
  @Input() placeholderText: string = 'FASTQ or FASTA file';
  @Input() showUploadNote: boolean = true;
  @Input() uploadNotePath: string = '#';
  @Input() uploadSpeed: number = 20; // milliseconds per 10%
  @Input() buttonClass: string = 'browse-button';
  
  @Output() fileSelected = new EventEmitter<FileUploadResult>();
  @Output() uploadComplete = new EventEmitter<FileUploadResult>();
  @Output() uploadProgress = new EventEmitter<number>();

  private apiService = inject(ApiService);

  selectedFile: File | null = null;
  fileName: string = '';
  savedFileName: string = '';
  isComplete: boolean = false;
  progress: number = 0;
  isUploading: boolean = false;
  uploadError: string = '';
  
  // Generate unique ID for each upload component instance
  readonly uploadId: string = `fileUpload-${Math.random().toString(36).substr(2, 9)}`;

  ngOnInit(): void {
    this.fileName = this.placeholderText;
  }

  onFileSelected(event: any): void {
    const file = event.target.files[0];
    if (file) {
      this.selectedFile = file;
      this.fileName = file.name;
      this.isComplete = false;
      this.progress = 0;
      this.isUploading = true;
      this.uploadError = '';

      // Emit file selected event
      this.fileSelected.emit({
        file: file,
        fileName: file.name,
        uploadComplete: false
      });

      // Start visual progress simulation
      const progressInterval = setInterval(() => {
        if (this.progress < 90) {
          this.progress += 10;
          this.uploadProgress.emit(this.progress);
        }
      }, this.uploadSpeed);

      // Upload file to backend
      this.apiService.uploadFile(file).subscribe({
        next: (response) => {
          clearInterval(progressInterval);
          this.progress = 100;
          this.isComplete = true;
          this.isUploading = false;
          this.savedFileName = response.saved_filename;
          
          // Emit upload complete event
          this.uploadComplete.emit({
            file: file,
            fileName: file.name,
            savedFileName: response.saved_filename,
            uploadComplete: true
          });
        },
        error: (error) => {
          clearInterval(progressInterval);
          this.isUploading = false;
          this.progress = 0;
          this.uploadError = error.error?.detail || 'Upload failed';
          
          this.uploadComplete.emit({
            file: file,
            fileName: file.name,
            uploadComplete: false,
            error: this.uploadError
          });
        }
      });
    }
  }

  triggerFileInput(): void {
    const fileInput = document.getElementById(this.uploadId) as HTMLInputElement;
    fileInput?.click();
  }

  resetUpload(): void {
    this.selectedFile = null;
    this.fileName = this.placeholderText;
    this.savedFileName = '';
    this.isComplete = false;
    this.progress = 0;
    this.isUploading = false;
    this.uploadError = '';
  }

  getUploadResult(): FileUploadResult | null {
    if (!this.selectedFile) {
      return null;
    }
    return {
      file: this.selectedFile,
      fileName: this.fileName,
      savedFileName: this.savedFileName,
      uploadComplete: this.isComplete,
      error: this.uploadError
    };
  }
}
