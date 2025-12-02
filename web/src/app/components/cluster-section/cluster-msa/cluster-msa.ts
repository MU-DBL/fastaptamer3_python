import { ChangeDetectorRef, Component, inject, OnInit, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Table, TableConfig } from '../../common/table/table';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { PlotModal } from '../../common/plot-modal/plot-modal';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { catchError, finalize, of, switchMap, tap } from 'rxjs';
import { ColumnName, FileService } from '../../../shared/file-service';

@Component({
  selector: 'app-cluster-msa',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    PlotModal,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster-msa.html',
  styleUrl: './cluster-msa.scss',
})

export class ClusterMsa {

  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private plotModalService = inject(PlotModalService);

  private cdr = inject(ChangeDetectorRef);

  selectedFile: File | null = null;
  savedFileName: string = '';
  uploadComplete: boolean = false;

  isProcessing = signal(false);
  isProcessingEntropy = signal(false);
  isProcessingMutInfo = signal(false);
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
  availablePalettes = ['magma', 'viridis', 'plasma', 'inferno', 'cividis', 'rocket', 'mako', 'turbo'];

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

  onEntropyPlot(): void {
    if (!this.processedFileName()) {
      alert('Please run MSA first before generating entropy plot.');
      return;
    }

    this.isProcessingEntropy.set(true);
    
    const params = {
      input_path: this.processedFileName()
    };

    this.apiService.clusterMsaEntropy(params).pipe(
      tap(response => {
        if (response.status === 'ok') {
          this.createEntropyPlot(response);
        }
      }),
      catchError(error => {
        console.error('Entropy calculation failed:', error);
        alert('Failed to calculate entropy. Please try again.');
        return of(null);
      }),
      finalize(() => {
        this.isProcessingEntropy.set(false);
      })
    ).subscribe();
  }

  onMIPlot(): void {
    if (!this.processedFileName()) {
      alert('Please run MSA first before generating mutual information plot.');
      return;
    }

    this.isProcessingMutInfo.set(true);
    
    const params = {
      input_path: this.processedFileName()
    };

    this.apiService.clusterMsaMutInfo(params).pipe(
      tap(response => {
        if (response.status === 'ok') {
          this.createMutualInfoPlot(response);
        }
      }),
      catchError(error => {
        console.error('Mutual information calculation failed:', error);
        alert('Failed to calculate mutual information. Please try again.');
        return of(null);
      }),
      finalize(() => {
        this.isProcessingMutInfo.set(false);
      })
    ).subscribe();
  }

  private createEntropyPlot(response: any): void {
    const plotData = [{
      x: response.positions,
      y: response.entropy_values,
      type: 'bar',
      name: 'Entropy',
      marker: {
        color: this.entropyBarFillColor,
        line: {
          color: this.entropyBarOutlineColor,
          width: 1
        }
      },
      hovertemplate: 'Position: %{x}<br>Entropy: %{y:.3f} nats<extra></extra>'
    }];

    const layout = {
      title: this.entropyPlotTitle,
      xaxis: {
        title: this.entropyXAxis,
        type: 'linear'
      },
      yaxis: {
        title: this.entropyYAxis
      },
      showlegend: false,
      hovermode: 'closest'
    };

    const config = {
      responsive: true,
      displayModeBar: true,
      displaylogo: false,
      toImageButtonOptions: {
        format: 'svg' as const,
        filename: 'msa_entropy'
      }
    };

    this.plotModalService.openPlot({
      data: plotData,
      layout: layout,
      config: config
    });
  }

  private createMutualInfoPlot(response: any): void {
    const colorscales: { [key: string]: string } = {
      'magma': 'Magma',
      'viridis': 'Viridis',
      'plasma': 'Plasma',
      'inferno': 'Inferno',
      'cividis': 'Cividis',
      'rocket': 'Magma', // Fallback to Magma since rocket might not be available in Plotly.js
      'mako': 'Blues',
      'turbo': 'Turbo'
    };

    const plotData = [{
      z: response.mutual_info_matrix,
      x: response.positions,
      y: response.positions,
      type: 'heatmap',
      colorscale: colorscales[this.miFillPalette] || 'Magma',
      showscale: true,
      colorbar: {
        title: this.miLegendTitle
      },
      hovertemplate: 'Position 1: %{x}<br>Position 2: %{y}<br>MI: %{z:.3f}<extra></extra>'
    }];

    const layout = {
      title: this.miPlotTitle,
      xaxis: {
        title: this.miXAxis,
        type: 'linear'
      },
      yaxis: {
        title: this.miYAxis,
        type: 'linear'
      },
      width: 600,
      height: 600
    };

    const config = {
      responsive: true,
      displayModeBar: true,
      displaylogo: false,
      toImageButtonOptions: {
        format: 'svg' as const,
        filename: 'msa_mutual_information'
      }
    };

    this.plotModalService.openPlot({
      data: plotData,
      layout: layout,
      config: config
    });
  }
  
}
