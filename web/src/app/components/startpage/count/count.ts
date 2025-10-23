import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';

@Component({
  selector: 'app-count',
  imports: [ 
    CommonModule, 
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS],
  templateUrl: './count.html',
  styleUrl: './count.scss'
})
export class Count {
  selectedFile: File | null = null;
  fileName: string = 'FASTQ or FASTA file';
  normalizeValue: string = '1e+06';
  returnReverseComplement: string = 'no';
  downloadFormat: string = 'fasta';
  uploadComplete: boolean = false;
  uploadProgress: number = 0;
  isUploading: boolean = false;

  onFileSelected(result: FileUploadResult): void {
    this.selectedFile = result.file;
    console.log('File selected:', result.fileName);
  }

  onUploadComplete(result: FileUploadResult): void {
    this.uploadComplete = true;
    console.log('Upload complete:', result.fileName);
  }

  onStart(): void {
    if (!this.uploadComplete) {
      alert('Please wait until the file upload is complete.');
      return;
    }
    console.log('Starting analysis with:', {
      file: this.selectedFile,
      normalizeValue: this.normalizeValue,
      returnReverseComplement: this.returnReverseComplement,
      downloadFormat: this.downloadFormat
    });
    // Add your start logic here
  }

  onDownload(): void {
    if (!this.uploadComplete) {
      alert('Please complete the analysis first.');
      return;
    }
    console.log('Downloading results...');
    // Add your download logic here
  }
}
