import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';

export interface FileUploadResult {
  file: File;
  fileName: string;
  uploadComplete: boolean;
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

  selectedFile: File | null = null;
  fileName: string = '';
  isComplete: boolean = false;
  progress: number = 0;
  isUploading: boolean = false;

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

      // Emit file selected event
      this.fileSelected.emit({
        file: file,
        fileName: file.name,
        uploadComplete: false
      });

      // Simulate upload process with progress
      const interval = setInterval(() => {
        this.progress += 10;
        
        if (this.progress >= 100) {
          clearInterval(interval);
          this.isComplete = true;
          this.isUploading = false;
          
          // Emit upload complete event
          this.uploadComplete.emit({
            file: file,
            fileName: file.name,
            uploadComplete: true
          });
        }
      }, this.uploadSpeed);
    }
  }

  triggerFileInput(): void {
    const fileInput = document.getElementById('fileUploadInput') as HTMLInputElement;
    fileInput?.click();
  }

  resetUpload(): void {
    this.selectedFile = null;
    this.fileName = this.placeholderText;
    this.isComplete = false;
    this.progress = 0;
    this.isUploading = false;
  }

  getUploadResult(): FileUploadResult | null {
    if (!this.selectedFile) {
      return null;
    }
    return {
      file: this.selectedFile,
      fileName: this.fileName,
      uploadComplete: this.isComplete
    };
  }
}
