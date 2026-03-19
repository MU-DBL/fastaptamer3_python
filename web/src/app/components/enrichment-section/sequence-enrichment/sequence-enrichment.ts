import { Component, inject, signal, ChangeDetectorRef, NgZone, OnDestroy } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { FileService } from '../../../shared/file-service';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { Table, TableConfig } from '../../common/table/table';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';

@Component({
  selector: 'app-sequence-enrichment',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './sequence-enrichment.html',
  styleUrl: './sequence-enrichment.scss',
  standalone: true
})
export class SequenceEnrichment implements OnDestroy {
  // Services
  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private plotModalService = inject(PlotModalService);
  private cdr = inject(ChangeDetectorRef);
  private ngZone = inject(NgZone);

  // Signals
  isProcessing = signal(false);
  processedFileName = signal('');

  // File handling - File 1
  selectedFile1: File | null = null;
  savedFileName1: string = '';
  uploadComplete1: boolean = false;

  // File handling - File 2
  selectedFile2: File | null = null;
  savedFileName2: string = '';
  uploadComplete2: boolean = false;

  // Enrichment parameters
  keepNA: string = 'no';
  downloadFormat: string = 'csv';

  // Data State
  enrichmentData: any[] = [];

  // Table configuration
  tableConfig: TableConfig = {
    columns: [
      { key: 'sequences', label: 'Sequences' },
      { key: 'ID.a', label: 'ID.a' },
      { key: 'Rank.a', label: 'Rank.a', exact_match: true },
      { key: 'Reads.a', label: 'Reads.a', exact_match: true },
      { key: 'RPU.a', label: 'RPU.a', exact_match: true },
      { key: 'ID.b', label: 'ID.b' },
      { key: 'Rank.b', label: 'Rank.b', exact_match: true },
      { key: 'Reads.b', label: 'Reads.b', exact_match: true },
      { key: 'RPU.b', label: 'RPU.b', exact_match: true },
      { key: 'Enrichment', label: 'Enrichment', exact_match: true },
      { key: 'log2E', label: 'log2E', exact_match: true }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  // Histogram plot customization
  adjustHistogram: string = 'no';
  histXAxis: string = 'log2(Enrichment)';
  histYAxis: string = 'Unique sequence count';
  histTitle: string = 'log2E histogram';
  histBarOutline: string = '#000000';
  histBarFill: string = '#87CEEB';

  // RPU Scatter plot customization
  adjustRPUScatter: string = 'no';
  rpuScatterXAxis: string = 'RPU.a';
  rpuScatterYAxis: string = 'RPU.b';
  rpuScatterTitle: string = 'RPU.a vs RPU.b';
  rpuScatterPointColor: string = '#87CEEB';

  // RA plot customization
  adjustRAPlot: string = 'no';
  raPlotXAxis: string = 'Average log2(RPU)';
  raPlotYAxis: string = 'Fold change';
  raPlotTitle: string = 'Enrichment RA plot';
  raPlotPointColor: string = '#87CEEB';

  // Box plot customization
  adjustBoxPlot: string = 'no';
  boxPlotXAxis: string = 'Cluster';
  boxPlotYAxis: string = 'Enrichment';
  boxPlotTitle: string = 'Enrichment distribution per cluster';
  boxPlotOutline: string = '#000000';
  boxPlotFill: string = '#87ceeb';

  // Watch for changes in adjustment toggles
  onAdjustHistogramChange(): void {
    if (this.adjustHistogram === 'no') {
      this.histXAxis = 'log2(Enrichment)';
      this.histYAxis = 'Unique sequence count';
      this.histTitle = 'log2E histogram';
      this.histBarOutline = '#000000';
      this.histBarFill = '#87CEEB';
    }
  }

  onAdjustRPUScatterChange(): void {
    if (this.adjustRPUScatter === 'no') {
      this.rpuScatterXAxis = 'RPU.a';
      this.rpuScatterYAxis = 'RPU.b';
      this.rpuScatterTitle = 'RPU.a vs RPU.b';
      this.rpuScatterPointColor = '#87CEEB';
    }
  }

  onAdjustRAPlotChange(): void {
    if (this.adjustRAPlot === 'no') {
      this.raPlotXAxis = 'Average log2(RPU)';
      this.raPlotYAxis = 'Fold change';
      this.raPlotTitle = 'Enrichment RA plot';
      this.raPlotPointColor = '#87CEEB';
    }
  }

  onAdjustBoxPlotChange(): void {
    if (this.adjustBoxPlot === 'no') {
      this.boxPlotXAxis = 'Cluster';
      this.boxPlotYAxis = 'Enrichment';
      this.boxPlotTitle = 'Enrichment distribution per cluster';
      this.boxPlotOutline = '#000000';
      this.boxPlotFill = '#87ceeb';
    }
  }

  // File 1 handlers
  onFile1Selected(result: FileUploadResult): void {
    if (this.savedFileName1) {
      this.apiService.deleteFile(this.savedFileName1).subscribe();
    }
    this.selectedFile1 = result.file;
    this.savedFileName1 = '';
    this.uploadComplete1 = false;
  }

  onFile1UploadComplete(result: FileUploadResult): void {
    this.savedFileName1 = result.savedFileName || '';
    this.uploadComplete1 = true;
  }

  // File 2 handlers
  onFile2Selected(result: FileUploadResult): void {
    if (this.savedFileName2) {
      this.apiService.deleteFile(this.savedFileName2).subscribe();
    }
    this.selectedFile2 = result.file;
    this.savedFileName2 = '';
    this.uploadComplete2 = false;
  }

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
  }

  ngOnDestroy(): void {
    if (this.savedFileName1) {
      this.apiService.deleteFile(this.savedFileName1).subscribe();
    }
    if (this.savedFileName2) {
      this.apiService.deleteFile(this.savedFileName2).subscribe();
    }
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
  }

  onFile2UploadComplete(result: FileUploadResult): void {
    this.savedFileName2 = result.savedFileName || '';
    this.uploadComplete2 = true;
  }

  onStart(): void {
    if (!this.uploadComplete1 || !this.uploadComplete2) {
      alert('Please upload both FASTA files before starting.');
      return;
    }

    if (!this.savedFileName1 || !this.savedFileName2) {
      alert('File upload incomplete. Please try again.');
      return;
    }

    this.isProcessing.set(true);
    this.enrichmentData = [];

    const params = {
      fadf1_cluster_path: this.savedFileName1,
      fadf2_cluster_path: this.savedFileName2,
      keep_na: this.keepNA === 'yes',
      output_format: this.downloadFormat
    };

    this.apiService.post('/sequence-enrich', params).pipe(
      tap(response => {
        console.log('Enrichment response:', response);
        this.processedFileName.set(response.result);
      }),
      switchMap(response => {
        if (response.result) {
          return this.apiService.downloadFile(response.result).pipe(
            tap(blob => this.parseFileBlob(blob, response.result))
          );
        }
        throw new Error('No result file returned from enrichment analysis');
      }),
      catchError(error => {
        const errorMessage = error.error?.detail || error.message || 'Unknown error';
        alert(`Enrichment failed: ${errorMessage}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
        this.cdr.detectChanges();
      })
    ).subscribe();
  }

  parseFileBlob(blob: Blob, filename: string): void {
    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      this.ngZone.run(() => {
        this.parseResultFile(content, filename);
      });
    };
    reader.readAsText(blob);
  }

  parseResultFile(content: string, filename: string): void {
    const isCsv = filename.endsWith('.csv');
    
    this.enrichmentData = [];
    
    if (isCsv) {
      // Parse CSV
      const lines = content.split('\n').filter(line => line.trim());
      if (lines.length === 0) return;
      
      const headers = lines[0].split(',').map(h => h.trim());
      
      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',').map(v => v.trim());
        const row: any = {};
        
        headers.forEach((header, index) => {
          const value = values[index];
          
          // Normalize column names - handle both 'sequences' and 'Sequences'
          let normalizedHeader = header;
          if (header === 'Sequences' || header === 'sequences') {
            normalizedHeader = 'sequences';
          }
          
          // Convert numeric fields
          if (['Rank.a', 'Rank.b', 'Reads.a', 'Reads.b', 'RPU.a', 'RPU.b', 'Enrichment', 'log2E', 'R', 'A'].includes(normalizedHeader)) {
            row[normalizedHeader] = value ? parseFloat(value) : 0;
          } else {
            row[normalizedHeader] = value || '';
          }
        });
        
        this.enrichmentData.push(row);
      }
    }
    
    console.log('Parsed enrichment data:', this.enrichmentData.length, 'rows');
    this.cdr.detectChanges();
  }

  onDownload(): void {
    if (!this.processedFileName()) {
      alert('No enrichment data available to download.');
      return;
    }

    this.apiService.downloadFile(this.processedFileName()).subscribe({
      next: (blob: Blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = this.processedFileName();
        link.click();
        window.URL.revokeObjectURL(url);
      },
      error: (error: any) => {
        console.error('Download error:', error);
        alert('Failed to download file.');
      }
    });
  }

  // ========================================================================
  // PLOTTING
  // ========================================================================

  async log2EnrichmentHistogram(): Promise<void> {
    if (this.enrichmentData.length === 0) {
      alert('No enrichment data available for plotting. Please run the analysis first.');
      return;
    }

    // Extract log2E values
    const log2EValues = this.enrichmentData
      .map(row => row['log2E'])
      .filter(val => val !== null && val !== undefined && isFinite(val));

    if (log2EValues.length === 0) {
      alert('No valid log2E values found for plotting.');
      return;
    }

    // Create histogram trace
    const trace = {
      x: log2EValues,
      type: 'histogram',
      nbinsx: 30,
      marker: {
        color: this.histBarFill,
        line: {
          color: this.histBarOutline,
          width: 1
        }
      },
      name: 'log2(Enrichment)'
    };

    const layout = {
      title: {
        text: this.histTitle,
        font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
      },
      xaxis: {
        title: {
          text: this.histXAxis,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0'
      },
      yaxis: {
        title: {
          text: this.histYAxis,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0'
      },
      autosize: true,
      height: 500,
      margin: { t: 60, b: 70, l: 90, r: 50 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white'
    };

    this.plotModalService.openPlot({
      data: [trace],
      layout: layout,
      config: { responsive: true }
    });
  }

  async rpuScatterPlot(): Promise<void> {
    if (this.enrichmentData.length === 0) {
      alert('No enrichment data available for plotting. Please run the analysis first.');
      return;
    }

    const EPSILON = 0.001;

    // Extract RPU values and sequences
    const rpuA = this.enrichmentData.map(row => (row['RPU.a'] || 0) + EPSILON);
    const rpuB = this.enrichmentData.map(row => (row['RPU.b'] || 0) + EPSILON);
    const sequences = this.enrichmentData.map(row => row['sequences'] || '');

    const trace = {
      x: rpuA,
      y: rpuB,
      mode: 'markers',
      type: 'scatter',
      marker: {
        color: this.rpuScatterPointColor,
        size: 6,
        opacity: 0.5
      },
      text: sequences,
      hovertemplate: '%{text}<br>RPU.a: %{x:.4f}<br>RPU.b: %{y:.4f}<extra></extra>'
    };

    const layout = {
      title: {
        text: this.rpuScatterTitle,
        font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
      },
      xaxis: {
        title: {
          text: this.rpuScatterXAxis,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        type: 'log',
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0'
      },
      yaxis: {
        title: {
          text: this.rpuScatterYAxis,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        type: 'log',
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0'
      },
      autosize: true,
      height: 500,
      margin: { t: 60, b: 70, l: 90, r: 50 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white',
      hovermode: 'closest'
    };

    this.plotModalService.openPlot({
      data: [trace],
      layout: layout,
      config: { responsive: true }
    });
  }

  async raPlot(): Promise<void> {
    if (this.enrichmentData.length === 0) {
      alert('No enrichment data available for plotting. Please run the analysis first.');
      return;
    }

    const EPSILON = 0.001;

    // Calculate R and A values
    const rValues: number[] = [];
    const aValues: number[] = [];
    const sequences: string[] = [];

    this.enrichmentData.forEach(row => {
      const rpuA = (row['RPU.a'] || 0) + EPSILON;
      const rpuB = (row['RPU.b'] || 0) + EPSILON;
      
      // R = log2(RPU.b / RPU.a)
      const r = Math.log2(rpuB / rpuA);
      // A = 0.5 * log2(RPU.b * RPU.a)
      const a = 0.5 * Math.log2(rpuB * rpuA);
      
      if (isFinite(r) && isFinite(a)) {
        rValues.push(r);
        aValues.push(a);
        sequences.push(row['sequences'] || '');
      }
    });

    if (rValues.length === 0) {
      alert('No valid data for RA plot.');
      return;
    }

    const trace = {
      x: aValues,
      y: rValues,
      mode: 'markers',
      type: 'scatter',
      marker: {
        color: this.raPlotPointColor,
        size: 6,
        opacity: 0.5
      },
      text: sequences,
      hovertemplate: '%{text}<br>A: %{x:.4f}<br>R: %{y:.4f}<extra></extra>'
    };

    const layout = {
      title: {
        text: this.raPlotTitle,
        font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
      },
      xaxis: {
        title: {
          text: this.raPlotXAxis,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0'
      },
      yaxis: {
        title: {
          text: this.raPlotYAxis,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0'
      },
      autosize: true,
      height: 500,
      margin: { t: 60, b: 70, l: 90, r: 50 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white',
      hovermode: 'closest'
    };

    this.plotModalService.openPlot({
      data: [trace],
      layout: layout,
      config: { responsive: true }
    });
  }

  async enrichmentBoxPlot(): Promise<void> {
    if (this.enrichmentData.length === 0) {
      alert('No enrichment data available for plotting. Please run the analysis first.');
      return;
    }

    // Check if cluster data exists
    const hasCluster = this.enrichmentData.some(row => row['Cluster.a'] !== undefined && row['Cluster.a'] !== null && row['Cluster.a'] !== 0 && row['Cluster.a'] !== '');

    let traces: any[];

    if (hasCluster) {
      // Group enrichment values by Cluster.a
      const clusterGroups = new Map<string, number[]>();
      this.enrichmentData.forEach(row => {
        const cluster = String(row['Cluster.a'] ?? 'N/A');
        const enrichment = parseFloat(row['Enrichment']);
        if (!isNaN(enrichment) && isFinite(enrichment)) {
          if (!clusterGroups.has(cluster)) clusterGroups.set(cluster, []);
          clusterGroups.get(cluster)!.push(enrichment);
        }
      });

      const sortedClusters = Array.from(clusterGroups.keys()).sort((a, b) => Number(a) - Number(b));
      traces = sortedClusters.map(cluster => ({
        type: 'box',
        y: clusterGroups.get(cluster),
        name: cluster,
        marker: { color: this.boxPlotFill, line: { color: this.boxPlotOutline, width: 1 } },
        boxmean: 'sd'
      }));
    } else {
      // Single box for all enrichment values
      const enrichmentValues = this.enrichmentData
        .map(row => parseFloat(row['Enrichment']))
        .filter(v => !isNaN(v) && isFinite(v));

      traces = [{
        type: 'box',
        y: enrichmentValues,
        name: 'All sequences',
        marker: { color: this.boxPlotFill, line: { color: this.boxPlotOutline, width: 1 } },
        boxmean: 'sd'
      }];
    }

    const layout = {
      title: { text: this.boxPlotTitle, font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' } },
      xaxis: { title: { text: this.boxPlotXAxis, font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' } }, showline: true, linewidth: 2, linecolor: 'black' },
      yaxis: { title: { text: this.boxPlotYAxis, font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' } }, showline: true, linewidth: 2, linecolor: 'black', showgrid: true, gridcolor: '#e0e0e0' },
      autosize: true,
      height: 500,
      margin: { t: 60, b: 70, l: 90, r: 50 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white',
      showlegend: false
    };

    this.plotModalService.openPlot({
      data: traces,
      layout: layout,
      config: { responsive: true }
    });
  }
}
