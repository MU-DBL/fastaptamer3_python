import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { Component, inject, signal } from '@angular/core';
import { ApiService, ProgressEvent } from '../../../shared/api.service';
import { FileService } from '../../../shared/file-service';

@Component({
  selector: 'app-cluster-phmm',
  imports: [  
    CommonModule, 
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './cluster-phmm.html',
  styleUrl: './cluster-phmm.scss',
})

export class ClusterPhmm {

  private apiService = inject(ApiService);
  private fileService = inject(FileService);

  // File upload state
  selectedFile: File | null = null;
  savedFileName: string = '';
  uploadComplete: boolean = false;
  
  // Use signals for reactive state
  isProcessing = signal<boolean>(false);

  // --- CHANGED: Separate signals for the two different output files ---
  phmmResultFile = signal<string>('');
  simulationResultFile = signal<string>('');

  // Slider parameters
  numSequences: number = 100;
  sequenceLength: number = 50;

  onFileSelected(result: FileUploadResult): void {
    this.selectedFile = result.file;
    
    // --- CHANGED: Reset both result signals when a new file is picked ---
    this.phmmResultFile.set('');
    this.simulationResultFile.set('');
    
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
    
    this.phmmResultFile.set('');
    this.simulationResultFile.set('');

    const params = {
      input_path: this.savedFileName,
      num_sequences: this.numSequences,
      sequence_length: this.sequenceLength,
      output_format_phmm: 'txt', 
      output_format_simulation: 'fasta',
      pseudocount_method: 'laplace'
    };

    this.apiService.clusterPhmmSimulate(params).subscribe({
      next: (response: any) => {
        console.log('Simulation successful:', response);
        this.isProcessing.set(false);

        if (response) {
          if (response.output_file_phmm) {
            this.phmmResultFile.set(response.output_file_phmm);
          }
          if (response.output_file_simulation) {
            this.simulationResultFile.set(response.output_file_simulation);
          }
        }
      },
      error: (error) => {
        console.error('Simulation failed:', error);
        this.isProcessing.set(false);
      }
    });
  }
  
  onDownloadPHMM(): void {
    const filename = this.phmmResultFile();
    if (!filename) {
      console.warn('No PHMM file available for download.');
      return;
    }
    // No need to replace/regex the name, the API gave us the exact name
    this.fileService.downloadFile(filename);
  }

  onDownloadSequences(): void {
    const filename = this.simulationResultFile();
    if (!filename) {
      console.warn('No simulation file available for download.');
      return;
    }
    this.fileService.downloadFile(filename);
  }
}