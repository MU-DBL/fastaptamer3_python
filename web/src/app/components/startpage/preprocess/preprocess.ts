import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { Component } from '@angular/core';

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
  selectedFile: File | null = null;
  fileName: string = 'FASTQ or FASTA file';
  downloadFormat: string = 'fasta';
  uploadComplete: boolean = false;
  isUploading: boolean = false;

  constant5Region: string = '';
  constant3Region: string = '';
  sequenceLengthMin: number = 10;
  sequenceLengthMax: number = 100;
  sequenceLengthRange: number = 500;
  maxAllowedError: number = 0.005;
  maxErrorRange: number = 1;

  onFileSelected(result: FileUploadResult): void {
    this.selectedFile = result.file;
    console.log('File selected:', result.fileName);
  }

  onUploadComplete(result: FileUploadResult): void {
    this.uploadComplete = true;
    console.log('Upload complete:', result.fileName);
  }

  onStart(): void {
    console.log('Starting with parameters:', {
      constant5Region: this.constant5Region,
      constant3Region: this.constant3Region,
      sequenceLength: [this.sequenceLengthMin, this.sequenceLengthMax],
      maxAllowedError: this.maxAllowedError
    });
  }

  onDownload(): void {
    console.log('Download triggered');
  }

  formatLabel(value: number): string {
    return `${value}`;
  }

  formatErrorLabel(value: number): string {
    return value.toFixed(3);
  }
}
