import { Component, inject, signal, ChangeDetectorRef, OnDestroy } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { SplitPanel } from '../../common/split-panel/split-panel';
import { catchError, finalize, switchMap, tap, map } from 'rxjs/operators';
import { of } from 'rxjs';
import { Table, TableConfig } from '../../common/table/table';
import { FileService } from '../../../shared/file-service';

@Component({
  selector: 'app-mutation-network',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    SplitPanel,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './mutation-network.html',
  styleUrl: './mutation-network.scss',
})


export class MutationNetwork implements OnDestroy {

  tableData: any[] = [];

  tableConfig: TableConfig = {
    columns: [
      { key: 'From_Sequence', label: 'From_Sequence' },
      { key: 'To_Sequence', label: 'To_Sequence' },
      { key: 'Transition_Cost', label: 'Transition_Cost' },
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  private cdr = inject(ChangeDetectorRef);
  private apiService = inject(ApiService);
  private fileService = inject(FileService);

  selectedFile: File | null = null;
  savedFileName: string = '';
  fileName: string = 'FASTA file';
  startSequence: string = '';
  endSequence: string = '';
  maxDistance: number = 1;
  uploadComplete: boolean = false;

  isProcessing = signal(false);
  processedFileName = signal('');
  errorMessage = signal('');

  get canStart(): boolean {
    return this.uploadComplete &&
      !this.isProcessing() &&
      this.startSequence.trim() !== '' &&
      this.endSequence.trim() !== '';
  }

  onFileSelected(result: FileUploadResult): void {
    if (this.savedFileName) {
      this.apiService.deleteFile(this.savedFileName).subscribe();
    }
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
    this.selectedFile = result.file;
    this.savedFileName = '';
    this.uploadComplete = false;
    this.processedFileName.set('');
    this.errorMessage.set('');
    this.tableData = [];
    console.log('File selected:', result.fileName);
  }

  ngOnDestroy(): void {
    if (this.savedFileName) {
      this.apiService.deleteFile(this.savedFileName).subscribe();
    }
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
  }

  onUploadComplete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete = true;
      this.savedFileName = result.savedFileName;
      console.log('Upload complete:', result.savedFileName);
    } else if (result.error) {
      console.error('Upload failed:', result.error);
      this.errorMessage.set('Upload failed: ' + result.error);
    }
  }

  onStart(): void {
    if (!this.uploadComplete || !this.savedFileName) {
      this.errorMessage.set('Please upload a file first!');
      return;
    }

    if (!this.startSequence || this.startSequence.trim() === '') {
      this.errorMessage.set('Please enter a start sequence!');
      return;
    }

    if (!this.endSequence || this.endSequence.trim() === '') {
      this.errorMessage.set('Please enter an end sequence!');
      return;
    }

    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.errorMessage.set('');
    this.tableData = [];

    const params = {
      input_path: this.savedFileName,
      start_node: this.startSequence.trim(),
      end_node: this.endSequence.trim(),
      max_cost: this.maxDistance,
      output_format: 'csv'
    };

    this.apiService.mutationNetwork(params).pipe(
      switchMap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Mutation Network completed:', response.result);

          return this.apiService.fetchFileText(response.result).pipe(
            map(text => this.fileService.parseResultFile(text, response.result)),
            tap(parsedData => {
              this.tableData = parsedData;
              this.cdr.detectChanges();
            })
          );
        }
        return of(null);
      }),
      catchError(error => {
        const errorMessage = error.error?.detail || error.message || 'Mutation network failed';
        alert(`Mutation network failed: ${errorMessage}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
  }

  onDownload(): void {
    if (!this.processedFileName()) {
      console.warn('No processed file available for download');
      return;
    }

    const link = document.createElement('a');
    link.href = this.apiService.getDownloadUrl(this.processedFileName());
    link.download = this.processedFileName();
    link.click();
  }
}
