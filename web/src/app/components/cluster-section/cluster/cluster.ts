import { Component, inject, signal, output, ChangeDetectorRef, OnDestroy } from '@angular/core';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { CommonModule } from '@angular/common';
import { SplitPanel } from '../../common/split-panel/split-panel';
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
    SplitPanel,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster.html',
  styleUrl: './cluster.scss',
})

export class Cluster implements OnDestroy {

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
  availableClusters: number[] = [];
  extractClusterNumber: number | null = null;

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
    this.clusterData = [];
    console.log('File selected:', result.fileName);
  }

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
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
            this.availableClusters = [...new Set(parsedData.map((r: any) => r[ColumnName.CLUSTER] as number))].sort((a, b) => a - b);
            this.extractClusterNumber = this.availableClusters[0] ?? null;
            this.cdr.detectChanges();
          })
        );
      }
      return of(null);
    }),
    catchError(error => {
        const errorMessage = error.error?.detail || error.message || 'Clustering failed';
        alert(`Clustering failed: ${errorMessage}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  onExtractCluster(): void {
    if (this.extractClusterNumber === null || this.clusterData.length === 0) return;

    const filtered = this.clusterData.filter(r => r[ColumnName.CLUSTER] === this.extractClusterNumber);
    if (filtered.length === 0) {
      alert(`No sequences found for cluster ${this.extractClusterNumber}.`);
      return;
    }

    const lines: string[] = [];
    filtered.forEach(r => {
      const header = `rank=${r[ColumnName.RANK]};read=${r[ColumnName.READS]};RPU=${r[ColumnName.RPU]};cluster=${r[ColumnName.CLUSTER]};RankInCluster=${r[ColumnName.RANK_IN_CLUSTER]};LED=${r[ColumnName.LED]}`;
      lines.push(`>${header}`);
      lines.push(r[ColumnName.SEQUENCES]);
    });

    const blob = new Blob([lines.join('\n')], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `cluster_${this.extractClusterNumber}.fasta`;
    a.click();
    URL.revokeObjectURL(url);
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
