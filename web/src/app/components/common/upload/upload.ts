import { Component, EventEmitter, Input, OnDestroy, Output, inject, ChangeDetectorRef, PLATFORM_ID } from '@angular/core';
import { CommonModule, isPlatformBrowser } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { Subscription } from 'rxjs';

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
export class Upload implements OnDestroy{
  @Input() acceptedFileTypes: string = '.fasta,.fastq,.fa,.fq';
  @Input() buttonText: string = 'Browse...';
  @Input() placeholderText: string = 'FASTQ or FASTA file';
  @Input() showUploadNote: boolean = true;
  @Input() uploadNotePath: string = '#';
  @Input() uploadSpeed: number = 200; // milliseconds per 10%
  @Input() buttonClass: string = 'browse-button';
  @Input() multiple: boolean = false;
  
  @Output() fileSelected = new EventEmitter<FileUploadResult>();
  @Output() uploadComplete = new EventEmitter<FileUploadResult>();
  @Output() uploadProgress = new EventEmitter<number>();

  private apiService = inject(ApiService);
  private cdr = inject(ChangeDetectorRef);
  private readonly platformId = inject(PLATFORM_ID);

  selectedFile: File | null = null;
  fileName: string = '';
  savedFileName: string = '';
  isComplete: boolean = false;
  progress: number = 0;
  isUploading: boolean = false;
  uploadError: string = '';

  // Multi-file mode state
  multiFiles: { name: string; progress: number; isComplete: boolean; error: string; savedFileName: string }[] = [];
  
  // Generate unique ID for each upload component instance
  readonly uploadId: string = `fileUpload-${Math.random().toString(36).substr(2, 9)}`;
  private progressInterval: any = null;
  private uploadSubscription: Subscription | null = null;
  private readonly onBeforeUnload = () => this.abortUpload();

  ngOnInit(): void {
    this.fileName = this.placeholderText;
    if (isPlatformBrowser(this.platformId)) {
      window.addEventListener('beforeunload', this.onBeforeUnload);
    }
  }

  ngOnDestroy(): void {
    if (isPlatformBrowser(this.platformId)) {
      window.removeEventListener('beforeunload', this.onBeforeUnload);
    }
    this.clearProgressInterval();
    this.abortUpload();
  }

  private abortUpload(): void {
    if (this.uploadSubscription) {
      this.uploadSubscription.unsubscribe();
      this.uploadSubscription = null;
    }
  }

  private clearProgressInterval(): void {
    if (this.progressInterval) {
      clearInterval(this.progressInterval);
      this.progressInterval = null;
    }
  }

  onFileSelected(event: any): void {
    if (this.multiple) {
      const files: File[] = Array.from(event.target.files || []);
      if (files.length === 0) return;

      // Reset previous uploads
      this.multiFiles = files.map(f => ({ name: f.name, progress: 0, isComplete: false, error: '', savedFileName: '' }));
      this.cdr.markForCheck();

      // Signal to parent that a new selection started (so it can clear its list)
      this.fileSelected.emit({ file: files[0], fileName: files[0].name, uploadComplete: false });

      files.forEach((file, i) => {
        const entry = this.multiFiles[i];
        const interval = setInterval(() => {
          if (entry.progress < 90) { entry.progress += 10; this.cdr.markForCheck(); }
        }, this.uploadSpeed);

        this.apiService.uploadFile(file).subscribe({
          next: (response) => {
            clearInterval(interval);
            entry.progress = 100;
            entry.isComplete = true;
            entry.savedFileName = response.saved_filename;
            this.cdr.markForCheck();
            this.uploadComplete.emit({ file, fileName: file.name, savedFileName: response.saved_filename, uploadComplete: true });
          },
          error: (error) => {
            clearInterval(interval);
            entry.progress = 0;
            entry.error = error.error?.detail || 'Upload failed';
            this.cdr.markForCheck();
            this.uploadComplete.emit({ file, fileName: file.name, uploadComplete: false, error: entry.error });
          }
        });
      });
      return;
    }

    const file = event.target.files[0];
    if (file) {
      this.clearProgressInterval();

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
      this.progressInterval = setInterval(() => {
        if (this.progress < 90) {
          this.progress += 1;
          this.uploadProgress.emit(this.progress);
          this.cdr.markForCheck();
        }
      }, this.uploadSpeed);

      // Upload file to backend
      this.uploadSubscription = this.apiService.uploadFile(file).subscribe({
        next: (response) => {
          this.clearProgressInterval();
          this.uploadSubscription = null;
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
          this.clearProgressInterval();
          this.uploadSubscription = null;
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
    this.clearProgressInterval();
    this.abortUpload();
    this.selectedFile = null;
    this.fileName = this.placeholderText;
    this.savedFileName = '';
    this.isComplete = false;
    this.progress = 0;
    this.isUploading = false;
    this.uploadError = '';
    const fileInput = document.getElementById(this.uploadId) as HTMLInputElement;
    if (fileInput) {
      fileInput.value = '';
    }
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
