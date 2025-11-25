import { Component, ViewChild, ElementRef, inject, signal, output, computed } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';
import { ChangeDetectorRef } from '@angular/core';

// Components & Services
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { ColumnName, FileService } from '../../../shared/file-service';
import { KmerAnalysisRequest, PlotlyTrace } from '../../../shared/kmer-analysis.types';
import { Table, TableConfig } from '../../common/table/table';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { ClusterColorService } from '../../../shared/cluster-color.service';

export interface DiversityResultsEvent {
  data: any[];
  inputFile?: string;
}

@Component({
  selector: 'app-cluster-diversity',
  standalone: true,
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './cluster-diversity.html',
  styleUrl: './cluster-diversity.scss',
})
export class ClusterDiversity {

  constructor(private cdr: ChangeDetectorRef, private plotModalService: PlotModalService) { }
  // Services
  private apiService = inject(ApiService);
  private fileService = inject(FileService);

  // DOM Elements (Replaces getElementById)
  @ViewChild('metaplotContainer') metaplotContainer!: ElementRef;
  @ViewChild('kmerPlotContainer') kmerPlotContainer!: ElementRef;

  // Outputs
  resultsReady = output<DiversityResultsEvent>();

  // Signals
  isProcessing = signal(false);
  isKmerProcessing = signal(false);
  processedFileName = signal('');

  // Plotting State
  plotly: any;

  metaplotData: any = null;
  metaplotParams: any = null;

  kmerPlotData: any = null;
  kmerPlotParams: any = null;

  // UI Flags
  showMetaplotModalFlag = false;
  showKmerPlotModalFlag = false;

  // Config - Metaplot
  adjustClusterMetadata = 'no';
  clusterMetaXAxis = "Cluster";
  clusterMetaPlotTitle = "Cluster metaplots";
  clusterMetaYAxis1 = "Seq. count";
  clusterMetaColor1 = "#0000ff";
  clusterMetaYAxis2 = "Read count";
  clusterMetaColor2 = "#ffa500";
  clusterMetaYAxis3 = "Avg. LED";
  clusterMetaColor3 = "#228b22";

  // Config - Kmer Plot
  kmerSize = "3";
  plotType = "UMAP";
  adjustKmerPlot = 'no';
  kmerPlotXAxis = "Dim1";
  kmerPlotYAxis = "Dim2";
  kmerPlotLegendTitle = "Cluster";
  kmerPlotTitle = "Cluster k-mer plot";
  kmerPlotColorPalette = "Dark2";

  // Data State
  selectedFile: File | null = null;
  savedFileName = '';
  uploadComplete = false;
  outputFormat = 'csv';

  diversityData: any[] = [];
  availableClusters = signal<number[]>([]); // All clusters from file
  selectedClusters = signal<number[]>([]);

  tableConfig: TableConfig = {
    columns: [
      { key: ColumnName.CLUSTER, label: 'Cluster', exact_match: true },
      { key: ColumnName.TOTAL_SEQUENCES, label: 'Total Sequences', exact_match: true },
      { key: ColumnName.TOTAL_READS, label: 'Total Reads', exact_match: true },
      { key: ColumnName.TOTAL_RPU, label: 'Total RPU', exact_match: true },
      { key: ColumnName.AVERAGE_LED, label: 'Average LED', exact_match: true },
      { key: ColumnName.SID, label: 'SeedID', exact_match: true }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  // ========================================================================
  // FILE HANDLING
  // ========================================================================

  onFileSelected(result: FileUploadResult): void {
    this.selectedFile = result.file;
    this.processedFileName.set('');
    this.diversityData = [];
    this.availableClusters = signal<number[]>([]); // All clusters from file
    this.selectedClusters = signal<number[]>([]);
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

  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) {
      console.warn('No file available for download.');
      return;
    }
    this.fileService.downloadFile(filename);
  }

  onStart(): void {
    if (!this.uploadComplete || !this.savedFileName) {
      console.warn('Please upload a clustered FASTA file first!');
      return;
    }

    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.diversityData = [];

    const params = {
      input_path: this.savedFileName,
      output_format: this.outputFormat
    };

    this.apiService.clusterDiversity(params).pipe(
      switchMap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);

          // Chain download and parsing
          return this.apiService.downloadFile(response.result).pipe(
            switchMap(blob => this.fileService.parseClusterFile(blob, response.result)),
            tap(parsedData => {
              this.diversityData = parsedData;
              const rawClusters = this.diversityData.map(row => row[ColumnName.CLUSTER]);
              const uniqueClusters = [...new Set(rawClusters)].sort((a, b) => a - b);
              this.availableClusters.set(uniqueClusters);
              console.log('Diversity data:', this.availableClusters());
            })
          );
        }
        return of(null);
      }),
      catchError(error => {
        alert('Clustering diversity failed!');
        return of(null);
      }),
      finalize(() => this.isProcessing.set(false))
    ).subscribe();
  }

  // ========================================================================
  // CLUSTER SELECTION HELPERS
  // ========================================================================

  addSelectedCluster(cluster: number): void {
    this.selectedClusters.update(current => {
      if (current.includes(cluster)) return current;
      return [...current, cluster].sort((a, b) => a - b);
    });
  }

  removeSelectedCluster(cluster: number): void {
    this.selectedClusters.update(current =>
      current.filter(c => c !== cluster)
    );
  }

  // If you need this for an event
  onClusterSelectionChange(clusters: number[]): void {
    this.selectedClusters.set(clusters);
  }

  // ========================================================================
  // METAPLOT LOGIC
  // ========================================================================

  clusterMetaPlot(): void {
    if (this.diversityData.length === 0) {
      alert('No data available. Run analysis first.');
      return;
    }

    const clusterNumbers = this.diversityData.map((d: any) => d[ColumnName.CLUSTER]);

    const plotConfig = [
      { key: ColumnName.TOTAL_SEQUENCES, axisY: 'y1', color: this.clusterMetaColor1, name: this.clusterMetaYAxis1 },
      { key: ColumnName.TOTAL_READS, axisY: 'y2', color: this.clusterMetaColor2, name: this.clusterMetaYAxis2 },
      { key: ColumnName.AVERAGE_LED, axisY: 'y3', color: this.clusterMetaColor3, name: this.clusterMetaYAxis3 }
    ];

    // Build Plotly traces
    const traces = plotConfig.map(cfg => ({
      x: clusterNumbers,
      y: this.diversityData.map((d: any) => d[cfg.key]),
      type: 'scatter' as const,
      mode: 'lines' as const,
      name: cfg.name,
      line: { color: cfg.color, width: 2 },
      xaxis: `x${cfg.axisY.replace('y', '')}`,  // x1, x2, x3
      yaxis: cfg.axisY
    }));

    // Build Plotly layout - CORRECTED!
    const layout = {
      title: {
        text: `${this.clusterMetaPlotTitle}`,
        font: { size: 16 },
        x: 0.5
      },
      showlegend: false,
      plot_bgcolor: 'white',

      // Top Plot (y1)
      xaxis: {
        domain: [0, 1],
        anchor: 'y1',
        showticklabels: false,
        tickfont: { size: 12 }
      },
      yaxis: {
        title: {
          text: this.clusterMetaYAxis1,
          font: { color: this.clusterMetaColor1 }
        },
        domain: [0.7, 1],  // Top 30% of plot
        tickfont: { color: this.clusterMetaColor1 }
      },

      // Middle Plot (y2)
      xaxis2: {
        domain: [0, 1],
        anchor: 'y2',
        showticklabels: false,
        tickfont: { size: 12 }
      },
      yaxis2: {
        title: {
          text: this.clusterMetaYAxis2,
          font: { color: this.clusterMetaColor2 }
        },
        domain: [0.35, 0.65],  // Middle 30% of plot
        tickfont: { color: this.clusterMetaColor2 }
      },

      // Bottom Plot (y3)
      xaxis3: {
        title: { text: this.clusterMetaXAxis },
        domain: [0, 1],
        anchor: 'y3',
        tickfont: { size: 12 }
      },
      yaxis3: {
        title: {
          text: this.clusterMetaYAxis3,
          font: { color: this.clusterMetaColor3 }
        },
        domain: [0, 0.3],  // Bottom 30% of plot
        tickfont: { color: this.clusterMetaColor3 }
      }
    };

    // Open modal with complete config
    this.plotModalService.openPlot({
      data: traces,
      layout: layout,
      config: {
        responsive: true,
        displayModeBar: true,
        toImageButtonOptions: {
          filename: 'cluster_metaplots',
          format: 'svg'
        }
      }
    });
  }
  // ========================================================================
  // K-MER PLOT LOGIC
  // ========================================================================
  kmerPlot(): void {
    if (this.diversityData.length === 0) {
      alert('No data available. Run analysis first.');
      return;
    }
    if (this.selectedClusters().length === 0) {
      alert('Please select at least one cluster.');
      return;
    }

    if (!this.savedFileName) {
      alert('No input file available.');
      return;
    }

    // Prepare request body
    const requestBody: KmerAnalysisRequest = {
      input_path: this.savedFileName,
      selected_clusters: this.selectedClusters(),
      k_size: parseInt(this.kmerSize),
      method: this.plotType.toLowerCase() as 'pca' | 'umap'
    };

    // Call API and render via modal service
    this.performKmerAnalysis(requestBody);
  }

  private performKmerAnalysis(requestBody: KmerAnalysisRequest): void {
    this.isKmerProcessing.set(true);

    this.apiService.clusterKmerAnalysis(requestBody).pipe(
      finalize(() => this.isKmerProcessing.set(false))
    ).subscribe({
      next: (response: any) => {
        if (response.status === 'ok' && response.data) {
          this.renderKmerPlotInModal(response.data, requestBody);
        } else {
          this.showKmerError('Invalid response from analysis');
        }
      },
      error: (error: any) => {
        this.showKmerError(`Analysis failed: ${error.error?.detail || error.message}`);
      }
    });
  }

  private renderKmerPlotInModal(data: any, requestBody: KmerAnalysisRequest): void {
    // 1. Group data by cluster (O(n) complexity)
    console.log('K-mer plot data received:', data);
    if (!data || !data.coordinates || !Array.isArray(data.coordinates)) {
      console.warn('K-mer analysis: No coordinates found in response.');
      return;
    }

    try {
      // Extract unique clusters and assign colors using color service
      const coordinates = data.coordinates;
      const clusters = coordinates.map((coord: any) => Number(coord.cluster));
      const uniqueSet = new Set(clusters);
      const uniqueClusters = Array.from(uniqueSet).sort((a, b) => (a as number) - (b as number)) as number[];
      const colorMap = ClusterColorService.createColorMap(uniqueClusters, this.kmerPlotColorPalette);

      // Group coordinates by cluster
      const clusterMap = new Map();
      data.coordinates.forEach((coord: any) => {
        const clusterId = Number(coord.cluster);

        if (!clusterMap.has(coord.cluster)) {
          clusterMap.set(clusterId, {
            x: [],
            y: [],
            sequences: [],
            color: colorMap[coord.cluster],
            cluster: coord.cluster
          });
        }
        const clusterData = clusterMap.get(clusterId);
        clusterData.x.push(coord.x);
        clusterData.y.push(coord.y);
        clusterData.sequences.push(coord.sequence || '');
      });

      // Create Plotly traces with proper typing
      const traces: PlotlyTrace[] = Array.from(clusterMap.values()).map(clusterData => ({
        x: clusterData.x,
        y: clusterData.y,
        mode: 'markers' as const,
        type: 'scatter' as const,
        name: `Cluster ${clusterData.cluster}`,
        marker: {
          color: clusterData.color,
          size: 8,
          opacity: 0.7
        },
        text: clusterData.sequences.map((seq: string) =>
          `Cluster: ${clusterData.cluster}<br>Sequence: ${seq.substring(0, 20)}...`
        ),
        textposition: 'top center' as const,
        hovertemplate: '<b>%{text}</b><br>X: %{x:.3f}<br>Y: %{y:.3f}<extra></extra>'
      }));

      // 3. Build layout
      const layout = {
        title: {
          text: `${this.kmerPlotTitle}`,
          font: { size: 16 },
          x: 0.5
        },
        xaxis: {
          title: this.kmerPlotXAxis,
          tickfont: { size: 12 }
        },
        yaxis: {
          title: this.kmerPlotYAxis,
          tickfont: { size: 12 }
        },
        legend: {
          title: { text: this.kmerPlotLegendTitle }
        },
        hovermode: 'closest' as const,
      };

      // 4. Open modal with plot config
      this.plotModalService.openPlot({
        data: traces,
        layout: layout,
        config: {
          responsive: true,
          displayModeBar: true,
          toImageButtonOptions: {
            filename: `kmer_${requestBody.method}_k${requestBody.k_size}`,
            format: 'svg'
          }
        }
      });
    }
   catch (error) {
        console.error("Error rendering K-mer plot:", error);
    }
}

  private showKmerError(message: string): void {
    this.isKmerProcessing.set(false);
    alert(message); // Or use a toast/snackbar service
  }
}