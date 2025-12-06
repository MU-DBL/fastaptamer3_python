import { Component, inject, signal, PLATFORM_ID, Inject, ChangeDetectorRef } from '@angular/core';
import { CommonModule, isPlatformBrowser } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';

// Components & Services
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { ColumnName, FileService } from '../../../shared/file-service';
import { Table, TableConfig } from '../../common/table/table';
import { PlotModalService } from '../../../shared/plot-modal.service';

@Component({
  selector: 'app-recluster',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './recluster.html',
  styleUrl: './recluster.scss',
})
export class Recluster {
  private cdr = inject(ChangeDetectorRef);
  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private plotModalService = inject(PlotModalService);

  constructor(@Inject(PLATFORM_ID) private platformId: Object) {}

  // ========================================================================
  // TABLE CONFIGURATION
  // ========================================================================
  tableConfig: TableConfig = {
    columns: [
      { key: ColumnName.SEQUENCES, label: 'Sequences' },
      { key: ColumnName.CLUSTER, label: 'Cluster', exact_match: true },
      { key: ColumnName.RANK_IN_CLUSTER, label: 'Rank In Cluster', exact_match: true },
      { key: ColumnName.LED, label: 'LED' },
      { key: ColumnName.ID_A, label: 'ID (Pop 1)' },
      { key: ColumnName.RANK_A, label: 'Rank (Pop 1)' },
      { key: ColumnName.READS_A, label: 'Reads (Pop 1)' },
      { key: ColumnName.RPU_A, label: 'RPU (Pop 1)' },
      { key: ColumnName.ID_B, label: 'ID (Pop 2)' },
      { key: ColumnName.RANK_B, label: 'Rank (Pop 2)' },
      { key: ColumnName.READS_B, label: 'Reads (Pop 2)' },
      { key: ColumnName.RPU_B, label: 'RPU (Pop 2)' },
      { key: ColumnName.ENRICHMENT, label: 'Enrichment' },
      { key: ColumnName.LOG2E, label: 'log2(Enrichment)' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  reclusterData: any[] = [];

  // ========================================================================
  // FILE UPLOAD STATE
  // ========================================================================
  selectedFile1: File | null = null;
  savedFileName1: string = '';
  uploadComplete1: boolean = false;

  selectedFile2: File | null = null;
  savedFileName2: string = '';
  uploadComplete2: boolean = false;

  // ========================================================================
  // PROCESSING STATE
  // ========================================================================
  isProcessing = signal(false);
  isHeatmapProcessing = signal(false);
  processedFileName = signal('');

  // ========================================================================
  // PARAMETERS
  // ========================================================================
  ledThreshold: number = 7;

  // ========================================================================
  // HEATMAP CUSTOMIZATION
  // ========================================================================
  adjustHeatmap: string = 'no';
  heatmapXAxis: string = 'Population 1 clusters';
  heatmapYAxis: string = 'Population 2 clusters';
  heatmapLegend: string = 'LED';
  heatmapTitle: string = 'LED between cluster seeds';
  heatmapPalette: string = 'Magma';

  // ========================================================================
  // POPULATION SIZE PLOT CUSTOMIZATION
  // ========================================================================
  adjustPopSize: string = 'no';
  popSizeXAxis: string = 'Cluster';
  popSizeYAxis: string = 'Sequence count';
  popSizeLegend: string = 'Population';
  popSizeTitle: string = 'Sequence count per cluster';
  popSizeColor1: string = '#1b9e77';
  popSizeColor2: string = '#d95f02';

  // ========================================================================
  // RPU PLOT CUSTOMIZATION
  // ========================================================================
  adjustRPU: string = 'no';
  rpuXAxis: string = 'Cluster';
  rpuYAxis: string = 'Avg. RPU';
  rpuLegend: string = 'Population';
  rpuTitle: string = 'Avg. RPU per cluster';
  rpuColor1: string = '#1b9e77';
  rpuColor2: string = '#d95f02';

  // ========================================================================
  // LED PLOT CUSTOMIZATION
  // ========================================================================
  adjustLED: string = 'no';
  ledXAxis: string = 'Cluster';
  ledYAxis: string = 'Avg. LED';
  ledTitle: string = 'Avg. LED per cluster';
  ledBarOutline: string = '#000000';
  ledBarFill: string = '#87ceeb';

  // ========================================================================
  // ENRICHMENT BOX PLOT CUSTOMIZATION
  // ========================================================================
  adjustEnrichment: string = 'no';
  enrichmentXAxis: string = 'Cluster';
  enrichmentTitle: string = 'Sequence enrichment per super-cluster';
  enrichmentBoxOutline: string = '#000000';
  enrichmentBoxFill: string = '#87ceeb';

  // ========================================================================
  // FILE UPLOAD HANDLERS - FILE 1
  // ========================================================================
  onFile1Selected(result: FileUploadResult): void {
    this.selectedFile1 = result.file;
    this.processedFileName.set('');
  }

  onUpload1Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete1 = true;
      this.savedFileName1 = result.savedFileName;
    } else if (result.error) {
      console.error('Upload 1 failed:', result.error);
    }
  }

  // ========================================================================
  // FILE UPLOAD HANDLERS - FILE 2
  // ========================================================================
  onFile2Selected(result: FileUploadResult): void {
    this.selectedFile2 = result.file;
    this.processedFileName.set('');
  }

  onUpload2Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete2 = true;
      this.savedFileName2 = result.savedFileName;
    } else if (result.error) {
      console.error('Upload 2 failed:', result.error);
    }
  }

  // ========================================================================
  // DOWNLOAD HANDLER
  // ========================================================================
  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) {
      console.warn('No file available for download.');
      return;
    }
    this.fileService.downloadFile(filename);
  }

  // ========================================================================
  // HEATMAP GENERATION
  // ========================================================================
  onGenerateHeatmap(): void {
    if (!this.uploadComplete1 || !this.savedFileName1) {
      console.warn('Please upload file 1 first!');
      return;
    }
    if (!this.uploadComplete2 || !this.savedFileName2) {
      console.warn('Please upload file 2 first!');
      return;
    }

    this.isHeatmapProcessing.set(true);

    const params = {
      fadf1_cluster_path: this.savedFileName1,
      fadf2_cluster_path: this.savedFileName2,
      led_threshold: this.ledThreshold,
      use_parallel: false
    };

    this.apiService.getReclusterLedMatrix(params)
      .pipe(
        tap((response) => {
          this.generateHeatmapPlot(response);
        }),
        catchError((error) => {
          console.error('LED matrix generation failed:', error);
          alert(`Heatmap generation failed: ${error.error?.detail || error.message}`);
          return of(null);
        }),
        finalize(() => {
          this.isHeatmapProcessing.set(false);
          this.cdr.detectChanges();
        })
      )
      .subscribe();
  }

  private generateHeatmapPlot(response: any): void {
    if (!isPlatformBrowser(this.platformId)) return;

    const { led_matrix, p1_cluster_ids, p2_cluster_ids } = response;

    const trace: any = {
      type: 'heatmap',
      z: led_matrix,
      x: p2_cluster_ids,
      y: p1_cluster_ids,
      colorscale: this.heatmapPalette,
      colorbar: {
        title: { text: this.heatmapLegend }
      },
      hoverongaps: false
    };

    const layout: any = {
      title: { text: this.heatmapTitle },
      xaxis: {
        title: { text: this.heatmapXAxis },
        side: 'bottom'
      },
      yaxis: {
        title: { text: this.heatmapYAxis }
      },
      autosize: true
    };

    const config = {
      responsive: true,
      displayModeBar: true,
      displaylogo: false
    };

    this.plotModalService.openPlot({
      data: [trace],
      layout: layout,
      config: config
    });
  }

  // ========================================================================
  // RECLUSTER
  // ========================================================================
  onRecluster(): void {
    if (!this.uploadComplete1 || !this.savedFileName1) {
      console.warn('Please upload file 1 first!');
      return;
    }
    if (!this.uploadComplete2 || !this.savedFileName2) {
      console.warn('Please upload file 2 first!');
      return;
    }

    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.reclusterData = [];

    const params = {
      fadf1_cluster_path: this.savedFileName1,
      fadf2_cluster_path: this.savedFileName2,
      led_threshold: this.ledThreshold,
      output_format: 'csv'
    };

    this.apiService.recluster(params)
      .pipe(
        switchMap((response) => {
          this.processedFileName.set(response.result);
          
          // Download and parse the result file
          return this.apiService.downloadFile(response.result).pipe(
            switchMap((blob) => this.fileService.parseClusterFile(blob, response.result))
          );
        }),
        tap((parsedData) => {
          this.reclusterData = parsedData;
        }),
        catchError((error) => {
          console.error('Reclustering failed:', error);
          alert(`Reclustering failed: ${error.error?.detail || error.message}`);
          return of([]);
        }),
        finalize(() => {
          this.isProcessing.set(false);
          this.cdr.detectChanges();
        })
      )
      .subscribe();
  }

  // ========================================================================
  // POPULATION SIZE PLOT
  // ========================================================================
  onGeneratePopSizePlot(): void {
    if (this.reclusterData.length === 0) {
      console.warn('Please generate recluster data first!');
      return;
    }

    // Group by cluster and count sequences per population
    const clusterGroups = new Map<number, { pop1: number; pop2: number }>();
    
    this.reclusterData.forEach(row => {
      const cluster = row[ColumnName.CLUSTER];
      if (!clusterGroups.has(cluster)) {
        clusterGroups.set(cluster, { pop1: 0, pop2: 0 });
      }
      
      const counts = clusterGroups.get(cluster)!;
      if (row[ColumnName.ID_A] !== undefined && row[ColumnName.ID_A] !== null && row[ColumnName.ID_A] !== '') {
        counts.pop1++;
      }
      if (row[ColumnName.ID_B] !== undefined && row[ColumnName.ID_B] !== null && row[ColumnName.ID_B] !== '') {
        counts.pop2++;
      }
    });

    const clusters = Array.from(clusterGroups.keys()).sort((a, b) => a - b);
    const pop1Counts = clusters.map(c => clusterGroups.get(c)!.pop1);
    const pop2Counts = clusters.map(c => clusterGroups.get(c)!.pop2);

    const trace1: any = {
      type: 'bar',
      name: 'Population 1',
      x: clusters,
      y: pop1Counts,
      marker: {
        color: this.popSizeColor1
      }
    };

    const trace2: any = {
      type: 'bar',
      name: 'Population 2',
      x: clusters,
      y: pop2Counts,
      marker: {
        color: this.popSizeColor2
      }
    };

    const layout: any = {
      title: { text: this.popSizeTitle },
      xaxis: { title: { text: this.popSizeXAxis } },
      yaxis: { title: { text: this.popSizeYAxis } },
      barmode: 'group',
      legend: { title: { text: this.popSizeLegend } },
      autosize: true
    };

    this.plotModalService.openPlot({
      data: [trace1, trace2],
      layout: layout,
      config: { responsive: true, displaylogo: false }
    });
  }

  // ========================================================================
  // AVERAGE RPU PLOT
  // ========================================================================
  onGenerateRPUPlot(): void {
    if (this.reclusterData.length === 0) {
      console.warn('Please generate recluster data first!');
      return;
    }

    // Group by cluster and calculate average RPU per population
    const clusterGroups = new Map<number, { 
      pop1Rpus: number[]; 
      pop2Rpus: number[] 
    }>();
    
    this.reclusterData.forEach(row => {
      const cluster = row[ColumnName.CLUSTER];
      if (!clusterGroups.has(cluster)) {
        clusterGroups.set(cluster, { pop1Rpus: [], pop2Rpus: [] });
      }
      
      const rpus = clusterGroups.get(cluster)!;
      
      const rpuA = row[ColumnName.RPU_A];
      if (rpuA !== undefined && rpuA !== null && rpuA !== '') {
        rpus.pop1Rpus.push(parseFloat(rpuA));
      }
      
      const rpuB = row[ColumnName.RPU_B];
      if (rpuB !== undefined && rpuB !== null && rpuB !== '') {
        rpus.pop2Rpus.push(parseFloat(rpuB));
      }
    });

    const clusters = Array.from(clusterGroups.keys()).sort((a, b) => a - b);
    const avgRpuPop1 = clusters.map(c => {
      const rpus = clusterGroups.get(c)!.pop1Rpus;
      return rpus.length > 0 ? rpus.reduce((a, b) => a + b, 0) / rpus.length : 0;
    });
    const avgRpuPop2 = clusters.map(c => {
      const rpus = clusterGroups.get(c)!.pop2Rpus;
      return rpus.length > 0 ? rpus.reduce((a, b) => a + b, 0) / rpus.length : 0;
    });

    const trace1: any = {
      type: 'bar',
      name: 'Population 1',
      x: clusters,
      y: avgRpuPop1,
      marker: {
        color: this.rpuColor1
      }
    };

    const trace2: any = {
      type: 'bar',
      name: 'Population 2',
      x: clusters,
      y: avgRpuPop2,
      marker: {
        color: this.rpuColor2
      }
    };

    const layout: any = {
      title: { text: this.rpuTitle },
      xaxis: { title: { text: this.rpuXAxis } },
      yaxis: { title: { text: this.rpuYAxis } },
      barmode: 'group',
      legend: { title: { text: this.rpuLegend } },
      autosize: true
    };

    this.plotModalService.openPlot({
      data: [trace1, trace2],
      layout: layout,
      config: { responsive: true, displaylogo: false }
    });
  }

  // ========================================================================
  // AVERAGE LED PLOT
  // ========================================================================
  onGenerateLEDPlot(): void {
    if (this.reclusterData.length === 0) {
      console.warn('Please generate recluster data first!');
      return;
    }

    // Group by cluster and calculate average LED
    const clusterGroups = new Map<number, number[]>();
    
    this.reclusterData.forEach(row => {
      const cluster = row[ColumnName.CLUSTER];
      const led = row[ColumnName.LED];
      
      if (led !== undefined && led !== null && led !== '') {
        if (!clusterGroups.has(cluster)) {
          clusterGroups.set(cluster, []);
        }
        clusterGroups.get(cluster)!.push(parseFloat(led));
      }
    });

    const clusters = Array.from(clusterGroups.keys()).sort((a, b) => a - b);
    const avgLed = clusters.map(c => {
      const leds = clusterGroups.get(c)!;
      return leds.reduce((a, b) => a + b, 0) / leds.length;
    });

    const trace: any = {
      type: 'bar',
      x: clusters,
      y: avgLed,
      marker: {
        color: this.ledBarFill,
        line: {
          color: this.ledBarOutline,
          width: 1
        }
      }
    };

    const layout: any = {
      title: { text: this.ledTitle },
      xaxis: { title: { text: this.ledXAxis } },
      yaxis: { title: { text: this.ledYAxis } },
      autosize: true
    };

    this.plotModalService.openPlot({
      data: [trace],
      layout: layout,
      config: { responsive: true, displaylogo: false }
    });
  }

  // ========================================================================
  // ENRICHMENT BOX PLOT
  // ========================================================================
  onGenerateEnrichmentPlot(): void {
    if (this.reclusterData.length === 0) {
      console.warn('Please generate recluster data first!');
      return;
    }

    // Group enrichment values by cluster
    const clusterGroups = new Map<number, number[]>();
    
    this.reclusterData.forEach(row => {
      const cluster = row[ColumnName.CLUSTER];
      const enrichment = row[ColumnName.ENRICHMENT];
      
      if (enrichment !== undefined && enrichment !== null && enrichment !== '' && !isNaN(enrichment) && isFinite(enrichment)) {
        if (!clusterGroups.has(cluster)) {
          clusterGroups.set(cluster, []);
        }
        clusterGroups.get(cluster)!.push(parseFloat(enrichment));
      }
    });

    const clusters = Array.from(clusterGroups.keys()).sort((a, b) => a - b);
    const traces = clusters.map(cluster => ({
      type: 'box',
      y: clusterGroups.get(cluster),
      name: `${cluster}`,
      marker: {
        color: this.enrichmentBoxFill,
        line: {
          color: this.enrichmentBoxOutline,
          width: 1
        }
      },
      boxmean: 'sd'
    }));

    const layout: any = {
      title: { text: this.enrichmentTitle },
      xaxis: { title: { text: this.enrichmentXAxis } },
      yaxis: { title: { text: 'Enrichment' } },
      autosize: true,
      showlegend: false
    };

    this.plotModalService.openPlot({
      data: traces,
      layout: layout,
      config: { responsive: true, displaylogo: false }
    });
  }

}
