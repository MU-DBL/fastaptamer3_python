import { Component, inject, signal, output, ChangeDetectorRef } from '@angular/core';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';
import { Table, TableConfig } from '../../common/table/table';
import { ColumnName, FileService } from '../../../shared/file-service';

export interface ClusterResultsEvent {
  data: any[];
  inputFile?: string;
}

@Component({
  selector: 'app-cluster',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster.html',
  styleUrl: './cluster.scss',
})

export class Cluster {

  private cdr = inject(ChangeDetectorRef);
  private apiService = inject(ApiService);
  private fileService = inject(FileService);

  tableConfig: TableConfig = {
    columns: [
      { key: ColumnName.ID, label: 'ID'},
      { key: ColumnName.CLUSTER, label: 'Cluster', exact_match: true },
      { key: ColumnName.RANK_IN_CLUSTER, label: 'Rank In Cluster', exact_match: true },
      { key: ColumnName.LED, label: 'LED' },
      { key: ColumnName.READS, label: 'Reads' },
      { key: ColumnName.RANK, label: 'Rank' },
      { key: ColumnName.RPU, label: 'RPU' },
      { key: ColumnName.SEQUENCES, label: 'Sequence' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  clusterData: any[] = [];

  resultsReady = output<ClusterResultsEvent>();

  downloadFormat: string = 'fasta';
  keepNonClusteredSequence: string = 'no';

  isProcessing = signal(false);
  processedFileName = signal('');

  selectedFile: File | null = null;
  savedFileName: string = '';
  uploadComplete: boolean = false;

  number_reads_to_cluster: number = 10;
  led: number = 7;
  number_of_cluster: number = 20;

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
    this.clusterData = []; // Clear data at start

    const params = {
      input_path: this.savedFileName,
      output_format: this.downloadFormat,
      min_reads: this.number_reads_to_cluster,
      max_led: this.led,
      total_clusters: this.number_of_cluster,
      keep_nc: this.keepNonClusteredSequence === 'yes'
    };

    console.log('Starting clustering with params:', params);

    this.apiService.clusterLed(params).pipe(
    switchMap(response => {
      if (response.status === 'ok' && response.result) {
        this.processedFileName.set(response.result);
        console.log('Clustering completed:', response.result);
        
        return this.apiService.downloadFile(response.result).pipe(
          switchMap(blob => 
            this.fileService.parseClusterFile(blob, response.result)
          ),
          tap(parsedData => {
            this.clusterData = parsedData;
            this.cdr.detectChanges();
          })
        );
      }
      return of(null);
    }),
    catchError(error => {
        console.error('Clustering failed:', error);
        let errorMessage = 'Clustering failed!';
        alert(errorMessage);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) {
      console.warn('No file available for download. Please run clustering first.');
      return;
    }

    console.log('Downloading file:', filename);
    this.fileService.downloadFile(filename);
  }
}
