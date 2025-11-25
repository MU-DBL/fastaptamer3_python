import { ChangeDetectorRef, Component, inject, OnInit, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Table, TableConfig } from '../../common/table/table';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { catchError, finalize, of, switchMap, tap } from 'rxjs';
import { ColumnName, FileService } from '../../../shared/file-service';

@Component({
  selector: 'app-cluster-msa',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster-msa.html',
  styleUrl: './cluster-msa.scss',
})

export class ClusterMsa {

  private apiService = inject(ApiService);
  private fileService = inject(FileService);

  private cdr = inject(ChangeDetectorRef);

  selectedFile: File | null = null;
  savedFileName: string = '';
  uploadComplete: boolean = false;

  isProcessing = signal(false);
  processedFileName = signal('');

  availableClusters: number[] = []; // Populates the mat-select
  selectedCluster: number | null = null;
  sequenceType: 'nucleotide' | 'aminoacid' = 'nucleotide';
  msaDownloadFormat: 'fasta' | 'csv' = 'fasta';
  

  // Entropy plot parameters
  adjustEntropyPlot = 'no';
  entropyXAxis = 'MSA position';
  entropyYAxis = 'Entropy (nats)';
  entropyLegendTitle = 'Nats';
  entropyPlotTitle = 'Entropy by position in MSA';
  entropyBarOutlineColor = '#000000';
  entropyBarFillColor = '#87CEEB';

  // Mutual information plot parameters
  adjustMIPlot = 'no';
  miXAxis = 'MSA position';
  miYAxis = 'MSA position';
  miLegendTitle = 'MI';
  miPlotTitle = 'Pairwise mutual information in MSA';
  miFillPalette = 'magma';
  availablePalettes = ['magma', 'viridis', 'plasma', 'inferno', 'cividis'];

  processedFile = '';

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

  onFileSelected(result: FileUploadResult): void {
    this.selectedFile = result.file;
    this.processedFileName.set('');
    this.availableClusters = [];
    this.selectedCluster = null;
    console.log('File selected:', result.fileName);
  }

  onUploadComplete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete = true;
      this.savedFileName = result.savedFileName;
      console.log('Upload complete:', result.savedFileName);

      this.loadClusterList(); 
    } else if (result.error) {
      console.error('Upload failed:', result.error);
    }
  }

  loadClusterList(): void {
    this.isProcessing.set(true);
    this.apiService.getClusterList({ input_path: this.savedFileName }).subscribe({
      next: (response) => {
        if (response.status === 'ok' || response.clusters) {
          this.availableClusters = response.clusters;
          console.log('Available clusters loaded:', this.availableClusters);
          // Optional: Auto-select the first cluster
          if (this.availableClusters.length > 0) {
            this.selectedCluster = this.availableClusters[0];
          }
        }
        this.isProcessing.set(false);
      },
      error: (error) => {
        console.error('Failed to load clusters:', error);
        this.isProcessing.set(false);
        alert('File uploaded, but failed to extract cluster list.');
      }
    });
  }

  onMSAStart(): void {
    if (!this.uploadComplete || !this.savedFileName) {
      console.warn('Please upload a file first!');
      return;
    }

    if (this.selectedCluster === null) {
      alert('Please select a cluster first.');
      return;
    }

    this.isProcessing.set(true);
    this.processedFileName.set('');

    // Map the HTML values to API parameters
    const params = {
      input_path: this.savedFileName,
      output_format: this.msaDownloadFormat,
      seq_type: this.sequenceType === 'nucleotide' ? 'dna' : 'protein', // Mapping 'nucleotide' -> 'dna' based on your previous API definition
      cluster_selected: this.selectedCluster
    };

    console.log('Starting MSA with params:', params);

    this.apiService.clusterMsa(params).pipe(
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
        let errorMessage = 'Clustering MSA failed!';
        alert(errorMessage);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  onMSADownload(): void {
    const filename = this.processedFileName();
    if (!filename) {
      console.warn('No file available for download.');
      return;
    }

    this.fileService.downloadFile(filename)
  }

  onMIPlot(): void {
    const params = {
      cluster: this.selectedCluster,
      adjustPlot: this.adjustMIPlot === 'yes',
      xAxis: this.miXAxis,
      yAxis: this.miYAxis,
      legendTitle: this.miLegendTitle,
      plotTitle: this.miPlotTitle,
      fillPalette: this.miFillPalette
    };
  }

    onEntropyPlot(): void {
    const params = {
      cluster: this.selectedCluster,
      adjustPlot: this.adjustEntropyPlot === 'yes',
      xAxis: this.entropyXAxis,
      yAxis: this.entropyYAxis,
      legendTitle: this.entropyLegendTitle,
      plotTitle: this.entropyPlotTitle,
      barOutlineColor: this.entropyBarOutlineColor,
      barFillColor: this.entropyBarFillColor
    };
  }
  
}
