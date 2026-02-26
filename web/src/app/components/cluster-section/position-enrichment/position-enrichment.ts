import { ChangeDetectorRef, Component, inject, OnInit, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Table, TableConfig } from '../../common/table/table';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { catchError, finalize, of, switchMap, tap } from 'rxjs';
import { ColumnName, FileService } from '../../../shared/file-service';


@Component({
  selector: 'app-position-enrichment',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS],
  templateUrl: './position-enrichment.html',
  styleUrl: './position-enrichment.scss',
})
export class PositionEnrichment {
  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private plotModalService = inject(PlotModalService);
  private cdr = inject(ChangeDetectorRef);

  selectedFile: File | null = null;
  savedFileName: string = '';
  uploadComplete: boolean = false;

  isProcessing = signal(false);
  processedFileName = signal('');

  availableClusters: number[] = []; // Populates the mat-select
  selectedCluster: number | null = null;
  sequenceType: 'nucleotide' | 'aminoacid' = 'nucleotide';
  clusterData: any[] = [];

  adjustPositionEnrichmentPlot: 'yes' | 'no' = 'no';
  adjustHeatPlot: 'yes' | 'no' = 'no';

  // Position Enrichment Plot properties
  posEnrichXAxis: string = 'Position in Alignment';
  posEnrichYAxis: string = 'Mean Enrichment';
  posEnrichPlotTitle: string = 'Position Enrichment Analysis';
  posEnrichBarOutlineColor: string = '#000000';
  posEnrichBarFillColor: string = '#4CAF50';

  // Heatmap properties
  heatmapXAxis: string = 'Position in Alignment';
  heatmapYAxis: string = 'Possible Enrichment';
  heatmapLegendTitle: string = 'Mean Enrichment';
  heatmapPlotTitle: string = 'Enrichment Heatmap';
  heatPalletteColor: string = 'Magma';
  availablePalettes = ['Magma', 'Viridis', 'Plasma', 'Inferno', 'Cividis', 'Blues', 'Turbo'];
  positions: any;
  residues: any;
  enrichment_matrix: any;
  avg_enrichment: any;

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

  onStart(): void {
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
      fadf_recluster_path: this.savedFileName,
      output_format: 'csv',
      seq_type: this.sequenceType === 'nucleotide' ? 'dna' : 'protein', // Mapping 'nucleotide' -> 'dna' based on your previous API definition
      cluster_selected: this.selectedCluster
    };

    console.log('Starting with params:', params);

    this.apiService.getPositionEnrichment(params).pipe(
      switchMap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          this.positions = response.positions;
          this.residues = response.residues;
          this.avg_enrichment = response.avg_enrichment;
          this.enrichment_matrix = response.enrichment_matrix;
          console.log('Clustering completed:', response);
          return of(response);
        }
        return of(null);
      }),
      catchError(error => {
        let errorMessage = 'Position enrichment analysis failed!';
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
      console.warn('No file available for download.');
      return;
    }
    this.fileService.downloadFile(filename)
  }

  onPosEnrichPlot(): void {
    if (!this.processedFileName()) {
      console.warn('No processed file available');
      return;
    }

    const plotData = [{
      x: this.positions,
      y: this.avg_enrichment,
      type: 'bar',
      name: 'Position Enrichment',
      marker: {
        color: this.posEnrichBarFillColor,
        line: {
          color: this.posEnrichBarOutlineColor,
          width: 1
        }
      },
      hovertemplate: 'Position: %{x}<br>Avg. Enrichment: %{y:.3f}<extra></extra>'
    }];

    const layout = {
      title: this.posEnrichPlotTitle,
      xaxis: {
        title: this.posEnrichXAxis,
        type: 'linear'
      },
      yaxis: {
        title: this.heatmapYAxis
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
        filename: 'position_enrichment'
      }
    };

    this.plotModalService.openPlot({
      data: plotData,
      layout: layout,
      config: config
    });
  }

  onHeatPlot(): void {
    if (!this.processedFileName()) {
      console.warn('No processed file available');
      return;
    }
    const plotData = [{
      z: this.enrichment_matrix,
      x: this.positions,
      y: this.residues,
      type: 'heatmap',
      colorscale: this.heatPalletteColor || 'Magma',
      showscale: true,
      colorbar: {
        title: this.heatmapLegendTitle
      },
      hovertemplate: 'Position: %{x}<br>Residue: %{y}<br>Enrichment: %{z:.3f}<extra></extra>'
    }];

    const layout = {
      title: this.heatmapPlotTitle,
      xaxis: {
        title: this.heatmapXAxis,
        type: 'linear'
      },
      yaxis: {
        title: this.heatmapYAxis,
        type: 'category'
      }
    };

    const config = {
      responsive: true,
      displayModeBar: true,
      displaylogo: false,
      toImageButtonOptions: {
        format: 'svg' as const,
        filename: 'enrichment_heatmap'
      }
    };

    this.plotModalService.openPlot({
      data: plotData,
      layout: layout,
      config: config
    });
  }
}
