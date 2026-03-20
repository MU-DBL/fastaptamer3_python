import { ChangeDetectorRef, Component, inject, OnDestroy, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Table, TableConfig } from '../../common/table/table';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { SplitPanel } from '../../common/split-panel/split-panel';
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
    SplitPanel,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster-msa.html',
  styleUrl: './cluster-msa.scss',
})

export class ClusterMsa implements OnDestroy {

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
  availablePalettes = ['Magma', 'Viridis', 'Plasma', 'Inferno', 'Cividis', 'Blues', 'Turbo'];

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
    // Clean up previous session's files before starting a new one
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
    this.availableClusters = [];
    this.selectedCluster = null;
    this.clusterData = [];
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

  onLoadResult(): void {
    if (!this.savedFileName) return;
    this.clusterData = [];
    this.apiService.downloadFile(this.savedFileName).pipe(
      switchMap(blob => this.fileService.parseClusterFile(blob, this.savedFileName)),
      tap(parsedData => {
        this.clusterData = parsedData;
        this.processedFileName.set(this.savedFileName);
        this.cdr.detectChanges();
      })
    ).subscribe();
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
        const errorMessage = error.error?.detail || error.message || 'MSA failed';
        alert(`MSA failed: ${errorMessage}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  onAlignmentPlot(): void {
    if (!this.processedFileName() || this.clusterData.length === 0) {
      alert('Please run MSA first.');
      return;
    }
    this.createAlignmentGridPlot();
  }

  private createAlignmentGridPlot(): void {
    const rows = this.clusterData
      .map(row => ({
        id: (row[ColumnName.ID] as string) ?? 'Unknown',
        seq: (row[ColumnName.SEQUENCES] as string)?.toUpperCase() ?? ''
      }))
      .filter(r => r.seq.length > 0);

    if (rows.length === 0) return;

    const seqLen = rows[0].seq.length;
    const positions = Array.from({ length: seqLen }, (_, i) => i + 1);
    const rowNumbers = rows.map((_, i) => i + 1);

    const isProtein = this.sequenceType === 'aminoacid';

    let charToNum: Record<string, number>;
    let colorscale: [number, string][];
    let zmax: number;
    let hoverLabel: string;

    if (isProtein) {
      // 20 standard AAs ordered by chemical property group; index 0 = gap/unknown
      const aaOrder = ['A','V','L','I','M','F','W','Y','S','T','N','Q','K','R','H','D','E','C','G','P'];
      charToNum = { '-': 0 };
      aaOrder.forEach((aa, i) => { charToNum[aa] = i + 1; });

      // Colors grouped by chemical property:
      // gray=gap, orange=hydrophobic, purple=aromatic, green=polar, blue=positive, red=negative
      // gold=Cys, lightgray=Gly, pink=Pro
      const aaColors = [
        '#BDBDBD', // 0  gap / unknown
        '#FF9800', // 1  A  hydrophobic
        '#FF9800', // 2  V  hydrophobic
        '#FF9800', // 3  L  hydrophobic
        '#FF9800', // 4  I  hydrophobic
        '#FF9800', // 5  M  hydrophobic
        '#9C27B0', // 6  F  aromatic
        '#9C27B0', // 7  W  aromatic
        '#9C27B0', // 8  Y  aromatic
        '#4CAF50', // 9  S  polar
        '#4CAF50', // 10 T  polar
        '#4CAF50', // 11 N  polar
        '#4CAF50', // 12 Q  polar
        '#2196F3', // 13 K  positive
        '#2196F3', // 14 R  positive
        '#2196F3', // 15 H  positive
        '#F44336', // 16 D  negative
        '#F44336', // 17 E  negative
        '#FFD700', // 18 C  cysteine
        '#E0E0E0', // 19 G  glycine
        '#E91E63', // 20 P  proline
      ];

      // Build step colorscale: each integer i (0-20) gets its own discrete band
      zmax = aaColors.length - 1; // 20
      const step = 1 / zmax;
      colorscale = [];
      for (let i = 0; i < aaColors.length - 1; i++) {
        colorscale.push([i * step, aaColors[i]]);
        colorscale.push([(i + 1) * step - 0.001, aaColors[i]]);
      }
      colorscale.push([1.0, aaColors[aaColors.length - 1]]);

      hoverLabel = 'AA';
    } else {
      // Nucleotide: 0=gap, 1=A, 2=T, 3=G, 4=C
      charToNum = { 'A': 1, 'T': 2, 'G': 3, 'C': 4 };
      colorscale = [
        [0,     '#BDBDBD'], [0.249, '#BDBDBD'],
        [0.25,  '#4CAF50'], [0.499, '#4CAF50'],
        [0.5,   '#F44336'], [0.749, '#F44336'],
        [0.75,  '#FF9800'], [0.999, '#FF9800'],
        [1.0,   '#2196F3']
      ];
      zmax = 4;
      hoverLabel = 'Base';
    }

    // Build z (numeric) and text matrices — rows = sequences, cols = positions
    const zMatrix: number[][] = [];
    const textMatrix: string[][] = [];

    for (const row of rows) {
      const zRow: number[] = [];
      const textRow: string[] = [];
      for (let pos = 0; pos < seqLen; pos++) {
        const ch = row.seq[pos] ?? '-';
        zRow.push(charToNum[ch] ?? 0);
        textRow.push(ch === '-' ? '_' : ch);
      }
      zMatrix.push(zRow);
      textMatrix.push(textRow);
    }

    const heatmapTrace = {
      type: 'heatmap' as const,
      z: zMatrix,
      text: textMatrix,
      texttemplate: '%{text}',
      x: positions,
      y: rowNumbers,
      customdata: rows.map(r => Array(seqLen).fill(r.id)),
      colorscale,
      showscale: false,
      zmin: 0,
      zmax,
      hovertemplate: `ID: %{customdata}<br>Position: %{x}<br>${hoverLabel}: %{text}<extra></extra>`
    };

    const legendEntries = isProtein ? [
      { name: 'Gap / Unknown',        color: '#BDBDBD' },
      { name: 'Hydrophobic (A,V,L,I,M)', color: '#FF9800' },
      { name: 'Aromatic (F,W,Y)',     color: '#9C27B0' },
      { name: 'Polar (S,T,N,Q)',      color: '#4CAF50' },
      { name: 'Positive (K,R,H)',     color: '#2196F3' },
      { name: 'Negative (D,E)',       color: '#F44336' },
      { name: 'Cys (C)',              color: '#FFD700' },
      { name: 'Gly (G)',              color: '#E0E0E0' },
      { name: 'Pro (P)',              color: '#E91E63' },
    ] : [
      { name: 'Gap', color: '#BDBDBD' },
      { name: 'A',   color: '#4CAF50' },
      { name: 'T',   color: '#F44336' },
      { name: 'G',   color: '#FF9800' },
      { name: 'C',   color: '#2196F3' },
    ];

    const legendTraces = legendEntries.map(entry => ({
      type: 'scatter' as const,
      x: [null as any],
      y: [null as any],
      mode: 'markers' as const,
      name: entry.name,
      marker: { color: entry.color, size: 14, symbol: 'square' },
      showlegend: true,
    }));

    const plotData = [heatmapTrace, ...legendTraces];

    const cellHeight = Math.max(20, Math.min(40, 600 / rows.length));
    const plotHeight = rows.length * cellHeight + 120;

    const layout = {
      title: { text: 'MSA Alignment' },
      xaxis: { title: { text: 'MSA Position' }, automargin: true },
      yaxis: { title: { text: 'Sequence' }, automargin: true, tickfont: { size: 10 } },
      height: plotHeight,
      showlegend: true,
      legend: {
        title: { text: 'Color scheme' },
        orientation: 'v' as const,
        x: 1.02,
        y: 1,
        xanchor: 'left' as const,
      },
      margin: { r: 180 },
    };

    this.plotModalService.openPlot({
      data: plotData,
      layout,
      config: {
        responsive: true,
        displayModeBar: true,
        displaylogo: false,
        toImageButtonOptions: { format: 'svg' as const, filename: 'msa_alignment' }
      }
    });
  }

  ngOnDestroy(): void {
    if (this.savedFileName) {
      this.apiService.deleteFile(this.savedFileName).subscribe();
    }
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
  }

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
    this.isProcessingEntropy.set(false);
    this.isProcessingMutInfo.set(false);
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
        const errorMessage = error.error?.detail || error.message || 'Entropy calculation failed';
        alert(`Entropy failed: ${errorMessage}`);
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
        const errorMessage = error.error?.detail || error.message || 'Mutual information calculation failed';
        alert(`Mutual information failed: ${errorMessage}`);
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
      title: { text: this.entropyPlotTitle },
      xaxis: {
        title: { text: this.entropyXAxis },
        type: 'linear',
        automargin: true
      },
      yaxis: {
        title: { text: this.entropyYAxis },
        automargin: true
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
    const plotData = [{
      z: response.mutual_info_matrix,
      x: response.positions,
      y: response.positions,
      type: 'heatmap',
      colorscale: this.miFillPalette || 'Magma',
      showscale: true,
      colorbar: {
        title: { text: this.miLegendTitle }
      },
      hovertemplate: 'Position 1: %{x}<br>Position 2: %{y}<br>MI: %{z:.3f}<extra></extra>'
    }];

    const layout = {
      title: { text: this.miPlotTitle },
      xaxis: {
        title: { text: this.miXAxis },
        type: 'linear',
        automargin: true
      },
      yaxis: {
        title: { text: this.miYAxis },
        type: 'linear',
        automargin: true
      },
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
