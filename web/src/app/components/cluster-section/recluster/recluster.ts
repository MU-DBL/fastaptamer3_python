import { Component, inject, signal, PLATFORM_ID, Inject, ChangeDetectorRef, OnDestroy } from '@angular/core';
import { CommonModule, isPlatformBrowser } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { switchMap, tap, catchError, finalize, map } from 'rxjs/operators';
import { of, forkJoin } from 'rxjs';

import { FileUploadResult, Upload } from '../../common/upload/upload';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { SplitPanel } from '../../common/split-panel/split-panel';
import { ApiService } from '../../../shared/api.service';
import { ColumnName, FileService } from '../../../shared/file-service';
import { Table, TableConfig } from '../../common/table/table';
import { PlotModalService } from '../../../shared/plot-modal.service';

@Component({
  selector: 'app-recluster',
  imports: [CommonModule, FormsModule, Upload, Table, SplitPanel, ...MATERIAL_IMPORTS],
  templateUrl: './recluster.html',
  styleUrl: './recluster.scss',
})
export class Recluster implements OnDestroy {
  private cdr = inject(ChangeDetectorRef);
  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private plotModalService = inject(PlotModalService);

  constructor(@Inject(PLATFORM_ID) private platformId: Object) {}

  // ========================================================================
  // MODE
  // ========================================================================
  mode: 'two-pop' | 'multi-round' = 'two-pop';

  onModeChange(): void {
    this.reclusterData = [];
    this.processedFileName.set('');
    this.multiRoundLabels = [];
  }

  // ========================================================================
  // TABLE CONFIG (2-pop)
  // ========================================================================
  tableConfig: TableConfig = {
    columns: [
      { key: ColumnName.SUPER_CLUSTER, label: 'Super-cluster', exact_match: true },
      { key: ColumnName.SEED, label: 'Seed' },
      { key: ColumnName.SIZE_POP1, label: 'Size (Pop 1)' },
      { key: ColumnName.SIZE_POP2, label: 'Size (Pop 2)' },
      { key: ColumnName.AVG_RPU_POP1, label: 'Avg RPU (Pop 1)' },
      { key: ColumnName.AVG_RPU_POP2, label: 'Avg RPU (Pop 2)' },
      { key: ColumnName.ENRICHMENT, label: 'Enrichment (Avg)' },
      { key: ColumnName.LOG2E, label: 'log2E (Avg)' },
      { key: ColumnName.SEED_RPU_POP1, label: 'Seed RPU (Pop 1)' },
      { key: ColumnName.SEED_RPU_POP2, label: 'Seed RPU (Pop 2)' },
      { key: ColumnName.SEED_ENRICHMENT, label: 'Enrichment (Seed)' },
      { key: ColumnName.SEED_LOG2E, label: 'log2E (Seed)' },
      { key: ColumnName.STATUS, label: 'Status' },
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100],
  };

  // Dynamic table config for multi-round
  multiTableConfig: TableConfig = {
    columns: [],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100],
  };

  reclusterData: any[] = [];
  sequenceData: any[] = [];
  multiRoundLabels: string[] = [];

  sequenceTableConfig: TableConfig = {
    columns: [
      { key: ColumnName.SEQUENCES, label: 'Sequence' },
      { key: ColumnName.CLUSTER, label: 'Super-cluster', exact_match: true },
      { key: ColumnName.ORIGINAL_CLUSTER_A, label: 'Orig. Cluster (Pop 1)', exact_match: true },
      { key: ColumnName.ORIGINAL_CLUSTER_B, label: 'Orig. Cluster (Pop 2)', exact_match: true },
      { key: ColumnName.RANK_IN_CLUSTER, label: 'Rank In Cluster', exact_match: true },
      { key: ColumnName.LED, label: 'LED' },
      { key: ColumnName.RPU_A, label: 'RPU (Pop 1)' },
      { key: ColumnName.RPU_B, label: 'RPU (Pop 2)' },
      { key: ColumnName.ENRICHMENT, label: 'Enrichment' },
      { key: ColumnName.LOG2E, label: 'log2E' },
    ],
    initialPageSize: 25,
    pageSizeOptions: [10, 25, 50, 100],
  };

  // ========================================================================
  // FILE UPLOAD STATE - 2 POP
  // ========================================================================
  savedFileName1 = '';
  uploadComplete1 = false;
  savedFileName2 = '';
  uploadComplete2 = false;

  // ========================================================================
  // FILE UPLOAD STATE - MULTI ROUND
  // ========================================================================
  savedFileNameM1 = '';
  uploadCompleteM1 = false;
  savedFileNameM2 = '';
  uploadCompleteM2 = false;
  savedFileNameM3 = '';
  uploadCompleteM3 = false;
  round1Label = 'R1';
  round2Label = 'R2';
  round3Label = 'R3';

  // ========================================================================
  // PROCESSING STATE
  // ========================================================================
  isProcessing = signal(false);
  isHeatmapProcessing = signal(false);
  processedFileName = signal('');
  processedSequenceFileName = signal('');

  // ========================================================================
  // PARAMETERS
  // ========================================================================
  ledThreshold = 7;
  enrichmentType: 'avg' | 'seed' = 'avg';
  enrichmentThreshold = 1.0;

  // ========================================================================
  // HEATMAP CUSTOMIZATION
  // ========================================================================
  adjustHeatmap = 'no';
  heatmapXAxis = 'Population 1 clusters';
  heatmapYAxis = 'Population 2 clusters';
  heatmapLegend = 'LED';
  heatmapTitle = 'LED between cluster seeds';
  heatmapPalette = 'Magma';

  // ========================================================================
  // POPULATION SIZE PLOT CUSTOMIZATION
  // ========================================================================
  adjustPopSize = 'no';
  popSizeXAxis = 'Super-cluster';
  popSizeYAxis = 'Sequence count';
  popSizeLegend = 'Population';
  popSizeTitle = 'Sequence count per super-cluster';
  popSizeColor1 = '#1b9e77';
  popSizeColor2 = '#d95f02';

  // ========================================================================
  // RPU PLOT CUSTOMIZATION
  // ========================================================================
  adjustRPU = 'no';
  rpuXAxis = 'Super-cluster';
  rpuYAxis = 'Avg. RPU';
  rpuLegend = 'Population';
  rpuTitle = 'Avg. RPU per super-cluster';
  rpuColor1 = '#1b9e77';
  rpuColor2 = '#d95f02';

  // ========================================================================
  // ENRICHMENT PLOT CUSTOMIZATION
  // ========================================================================
  adjustEnrichment = 'no';
  enrichmentXAxis = 'Super-cluster';
  enrichmentYAxis = 'Enrichment';
  enrichmentTitle = 'Enrichment per super-cluster';
  enrichmentBarFill = '#87ceeb';
  enrichmentBarOutline = '#000000';

  // ========================================================================
  // RPU SCATTER PLOT CUSTOMIZATION (two-pop)
  // ========================================================================
  adjustRPUScatter = 'no';
  rpuScatterXAxis = 'Pop 1 Avg. RPU';
  rpuScatterYAxis = 'Pop 2 Avg. RPU';
  rpuScatterTitle = 'RPU scatter: Pop 1 vs Pop 2';
  rpuScatterColor = '#1f77b4';

  // ========================================================================
  // TRAJECTORY PLOT CUSTOMIZATION (multi-round)
  // ========================================================================
  adjustTrajectory = 'no';
  trajectoryYAxis = 'Avg. RPU';
  trajectoryTitle = 'Cluster trajectory across rounds';
  trajectoryTopN: number | null = null;

  // ========================================================================
  // LIFECYCLE
  // ========================================================================
  ngOnDestroy(): void {
    [this.savedFileName1, this.savedFileName2,
     this.savedFileNameM1, this.savedFileNameM2, this.savedFileNameM3,
     this.processedFileName(), this.processedSequenceFileName()].forEach(f => {
      if (f) this.apiService.deleteFile(f).subscribe();
    });
  }

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
  }

  // ========================================================================
  // FILE UPLOAD HANDLERS - 2 POP
  // ========================================================================
  onFile1Selected(result: FileUploadResult): void {
    if (this.savedFileName1) this.apiService.deleteFile(this.savedFileName1).subscribe();
    this.savedFileName1 = '';
    this.uploadComplete1 = false;
    this.processedFileName.set('');
  }

  onUpload1Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete1 = true;
      this.savedFileName1 = result.savedFileName;
    }
  }

  onFile2Selected(result: FileUploadResult): void {
    if (this.savedFileName2) this.apiService.deleteFile(this.savedFileName2).subscribe();
    this.savedFileName2 = '';
    this.uploadComplete2 = false;
    this.processedFileName.set('');
  }

  onUpload2Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete2 = true;
      this.savedFileName2 = result.savedFileName;
    }
  }

  // ========================================================================
  // FILE UPLOAD HANDLERS - MULTI ROUND
  // ========================================================================
  onFileM1Selected(result: FileUploadResult): void {
    if (this.savedFileNameM1) this.apiService.deleteFile(this.savedFileNameM1).subscribe();
    this.savedFileNameM1 = '';
    this.uploadCompleteM1 = false;
    this.processedFileName.set('');
  }

  onUploadM1Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadCompleteM1 = true;
      this.savedFileNameM1 = result.savedFileName;
    }
  }

  onFileM2Selected(result: FileUploadResult): void {
    if (this.savedFileNameM2) this.apiService.deleteFile(this.savedFileNameM2).subscribe();
    this.savedFileNameM2 = '';
    this.uploadCompleteM2 = false;
    this.processedFileName.set('');
  }

  onUploadM2Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadCompleteM2 = true;
      this.savedFileNameM2 = result.savedFileName;
    }
  }

  onFileM3Selected(result: FileUploadResult): void {
    if (this.savedFileNameM3) this.apiService.deleteFile(this.savedFileNameM3).subscribe();
    this.savedFileNameM3 = '';
    this.uploadCompleteM3 = false;
    this.processedFileName.set('');
  }

  onUploadM3Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadCompleteM3 = true;
      this.savedFileNameM3 = result.savedFileName;
    }
  }

  // ========================================================================
  // DOWNLOAD
  // ========================================================================
  onDownload(): void {
    const filename = this.processedFileName();
    if (filename) this.fileService.downloadFile(filename);
  }

  onDownloadSequences(): void {
    const filename = this.processedSequenceFileName();
    if (filename) this.fileService.downloadFile(filename);
  }

  // ========================================================================
  // HEATMAP (2-pop only)
  // ========================================================================
  onGenerateHeatmap(): void {
    if (!this.uploadComplete1 || !this.uploadComplete2) return;
    this.isHeatmapProcessing.set(true);

    this.apiService.getReclusterLedMatrix({
      fadf1_cluster_path: this.savedFileName1,
      fadf2_cluster_path: this.savedFileName2,
      led_threshold: this.ledThreshold,
      use_parallel: false,
    }).pipe(
      tap(response => this.generateHeatmapPlot(response)),
      catchError(error => {
        alert(`Heatmap generation failed: ${error.error?.detail || error.message}`);
        return of(null);
      }),
      finalize(() => { this.isHeatmapProcessing.set(false); this.cdr.detectChanges(); })
    ).subscribe();
  }

  private generateHeatmapPlot(response: any): void {
    if (!isPlatformBrowser(this.platformId)) return;
    const { led_matrix, p1_cluster_ids, p2_cluster_ids } = response;
    this.plotModalService.openPlot({
      data: [{
        type: 'heatmap',
        z: led_matrix,
        x: p2_cluster_ids,
        y: p1_cluster_ids,
        colorscale: this.heatmapPalette,
        colorbar: { title: { text: this.heatmapLegend } },
        hoverongaps: false,
      }],
      layout: {
        title: { text: this.heatmapTitle },
        xaxis: { title: { text: this.heatmapXAxis }, side: 'bottom' },
        yaxis: { title: { text: this.heatmapYAxis } },
        autosize: true,
      },
      config: { responsive: true, displayModeBar: true, displaylogo: false },
    });
  }

  // ========================================================================
  // RECLUSTER (2-pop)
  // ========================================================================
  onRecluster(): void {
    if (!this.uploadComplete1 || !this.uploadComplete2) return;
    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.processedSequenceFileName.set('');
    this.reclusterData = [];
    this.sequenceData = [];

    this.apiService.recluster({
      fadf1_cluster_path: this.savedFileName1,
      fadf2_cluster_path: this.savedFileName2,
      led_threshold: this.ledThreshold,
      enrichment_type: this.enrichmentType,
      enrichment_threshold: this.enrichmentThreshold,
      output_format: 'csv',
    }).pipe(
      switchMap((response: any) => {
        this.processedFileName.set(response.result);
        this.processedSequenceFileName.set(response.result_sequences ?? '');
        const clusters$ = this.apiService.fetchFileText(response.result).pipe(
          map(text => this.fileService.parseResultFile(text, response.result))
        );
        const sequences$ = response.result_sequences
          ? this.apiService.fetchFileText(response.result_sequences).pipe(
              map(text => this.fileService.parseResultFile(text, response.result_sequences))
            )
          : of([]);
        return forkJoin([clusters$, sequences$]);
      }),
      tap(([clusterData, seqData]) => {
        this.reclusterData = clusterData;
        this.sequenceData = seqData;
      }),
      catchError(error => {
        alert(`Reclustering failed: ${error.error?.detail || error.message}`);
        return of([[], []]);
      }),
      finalize(() => { this.isProcessing.set(false); this.cdr.detectChanges(); })
    ).subscribe();
  }

  // ========================================================================
  // RECLUSTER MULTI
  // ========================================================================
  onReclusterMulti(): void {
    if (!this.uploadCompleteM1 || !this.uploadCompleteM2 || !this.uploadCompleteM3) return;
    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.reclusterData = [];

    this.apiService.reclusterMulti({
      fadf1_cluster_path: this.savedFileNameM1,
      fadf2_cluster_path: this.savedFileNameM2,
      fadf3_cluster_path: this.savedFileNameM3,
      round1_label: this.round1Label,
      round2_label: this.round2Label,
      round3_label: this.round3Label,
      led_threshold: this.ledThreshold,
      output_format: 'csv',
    }).pipe(
      switchMap((response: any) => {
        this.processedFileName.set(response.result);
        this.multiRoundLabels = response.labels ?? [this.round1Label, this.round2Label, this.round3Label];
        this.buildMultiTableConfig();
        return this.apiService.fetchFileText(response.result).pipe(
          map(text => this.fileService.parseResultFile(text, response.result))
        );
      }),
      tap(parsedData => { this.reclusterData = parsedData; }),
      catchError(error => {
        alert(`Multi-round reclustering failed: ${error.error?.detail || error.message}`);
        return of([]);
      }),
      finalize(() => { this.isProcessing.set(false); this.cdr.detectChanges(); })
    ).subscribe();
  }

  private buildMultiTableConfig(): void {
    const [r1, r2, r3] = this.multiRoundLabels;
    this.multiTableConfig = {
      columns: [
        { key: 'SuperCluster', label: 'Super-cluster', exact_match: true },
        { key: 'Seed', label: 'Seed' },
        { key: `Size.${r1}`, label: `Size (${r1})` },
        { key: `Size.${r2}`, label: `Size (${r2})` },
        { key: `Size.${r3}`, label: `Size (${r3})` },
        { key: `AvgRPU.${r1}`, label: `Avg RPU (${r1})` },
        { key: `AvgRPU.${r2}`, label: `Avg RPU (${r2})` },
        { key: `AvgRPU.${r3}`, label: `Avg RPU (${r3})` },
        { key: `E.${r1}.${r2}`, label: `Enrichment (${r1}→${r2})` },
        { key: `E.${r2}.${r3}`, label: `Enrichment (${r2}→${r3})` },
        { key: `SeedRPU.${r1}`, label: `Seed RPU (${r1})` },
        { key: `SeedRPU.${r2}`, label: `Seed RPU (${r2})` },
        { key: `SeedRPU.${r3}`, label: `Seed RPU (${r3})` },
        { key: `SeedE.${r1}.${r2}`, label: `Seed Enrichment (${r1}→${r2})` },
        { key: `SeedE.${r2}.${r3}`, label: `Seed Enrichment (${r2}→${r3})` },
      ],
      initialPageSize: 10,
      pageSizeOptions: [10, 25, 50, 100],
    };
  }

  // ========================================================================
  // POPULATION SIZE PLOT
  // ========================================================================
  onGeneratePopSizePlot(): void {
    if (this.reclusterData.length === 0) return;
    const clusters = this.reclusterData.map(r => r[ColumnName.SUPER_CLUSTER]);
    const pop1 = this.reclusterData.map(r => r[ColumnName.SIZE_POP1] ?? 0);
    const pop2 = this.reclusterData.map(r => r[ColumnName.SIZE_POP2] ?? 0);

    this.plotModalService.openPlot({
      data: [
        { type: 'bar', name: 'Population 1', x: clusters, y: pop1, marker: { color: this.popSizeColor1 } },
        { type: 'bar', name: 'Population 2', x: clusters, y: pop2, marker: { color: this.popSizeColor2 } },
      ],
      layout: {
        title: { text: this.popSizeTitle },
        xaxis: { title: { text: this.popSizeXAxis } },
        yaxis: { title: { text: this.popSizeYAxis } },
        barmode: 'group',
        legend: { title: { text: this.popSizeLegend } },
        autosize: true,
      },
      config: { responsive: true, displaylogo: false },
    });
  }

  // ========================================================================
  // RPU PLOT
  // ========================================================================
  onGenerateRPUPlot(): void {
    if (this.reclusterData.length === 0) return;
    const clusters = this.reclusterData.map(r => r[ColumnName.SUPER_CLUSTER]);
    const rpu1 = this.reclusterData.map(r => r[ColumnName.AVG_RPU_POP1] ?? null);
    const rpu2 = this.reclusterData.map(r => r[ColumnName.AVG_RPU_POP2] ?? null);

    this.plotModalService.openPlot({
      data: [
        { type: 'bar', name: 'Population 1', x: clusters, y: rpu1, marker: { color: this.rpuColor1 } },
        { type: 'bar', name: 'Population 2', x: clusters, y: rpu2, marker: { color: this.rpuColor2 } },
      ],
      layout: {
        title: { text: this.rpuTitle },
        xaxis: { title: { text: this.rpuXAxis } },
        yaxis: { title: { text: this.rpuYAxis } },
        barmode: 'group',
        legend: { title: { text: this.rpuLegend } },
        autosize: true,
      },
      config: { responsive: true, displaylogo: false },
    });
  }

  // ========================================================================
  // ENRICHMENT BAR PLOT (2-pop)
  // ========================================================================
  onGenerateEnrichmentPlot(): void {
    if (this.reclusterData.length === 0) return;
    const clusters = this.reclusterData.map(r => r[ColumnName.SUPER_CLUSTER]);
    const enrichments = this.reclusterData.map(r => {
      const v = parseFloat(r[ColumnName.ENRICHMENT]);
      return isFinite(v) ? v : null;
    });

    this.plotModalService.openPlot({
      data: [{
        type: 'bar',
        x: clusters,
        y: enrichments,
        marker: { color: this.enrichmentBarFill, line: { color: this.enrichmentBarOutline, width: 1 } },
      }],
      layout: {
        title: { text: this.enrichmentTitle },
        xaxis: { title: { text: this.enrichmentXAxis } },
        yaxis: { title: { text: this.enrichmentYAxis } },
        autosize: true,
      },
      config: { responsive: true, displaylogo: false },
    });
  }

  // ========================================================================
  // RPU SCATTER PLOT (two-pop)
  // ========================================================================
  onGenerateRPUScatterPlot(): void {
    if (this.reclusterData.length === 0) return;
    const clusters = this.reclusterData.map(r => r[ColumnName.SUPER_CLUSTER]);
    const rpu1 = this.reclusterData.map(r => r[ColumnName.AVG_RPU_POP1] ?? null);
    const rpu2 = this.reclusterData.map(r => r[ColumnName.AVG_RPU_POP2] ?? null);
    const maxVal = Math.max(...rpu1.filter(v => v != null), ...rpu2.filter(v => v != null));

    this.plotModalService.openPlot({
      data: [
        {
          type: 'scatter',
          mode: 'markers+text',
          x: rpu1,
          y: rpu2,
          text: clusters.map((c: any) => `Cluster ${c}`),
          textposition: 'top center',
          marker: { color: this.rpuScatterColor, size: 8 },
        },
        {
          type: 'scatter',
          mode: 'lines',
          x: [0, maxVal],
          y: [0, maxVal],
          line: { color: '#888', dash: 'dash', width: 1 },
          showlegend: false,
          hoverinfo: 'none',
        },
      ],
      layout: {
        title: { text: this.rpuScatterTitle },
        xaxis: { title: { text: this.rpuScatterXAxis }, zeroline: false },
        yaxis: { title: { text: this.rpuScatterYAxis }, zeroline: false },
        autosize: true,
      },
      config: { responsive: true, displaylogo: false },
    });
  }

  // ========================================================================
  // TRAJECTORY LINE PLOT (multi-round)
  // ========================================================================
  onGenerateTrajectoryPlot(): void {
    if (this.reclusterData.length === 0 || this.multiRoundLabels.length === 0) return;
    const [r1, r2, r3] = this.multiRoundLabels;
    const rounds = [r1, r2, r3];

    const data = (this.trajectoryTopN && this.trajectoryTopN > 0)
      ? this.reclusterData.slice(0, this.trajectoryTopN)
      : this.reclusterData;

    const traces = data.map((row: any) => ({
      type: 'scatter',
      mode: 'lines+markers',
      name: `Cluster ${row['SuperCluster']}`,
      x: rounds,
      y: [row[`AvgRPU.${r1}`] ?? null, row[`AvgRPU.${r2}`] ?? null, row[`AvgRPU.${r3}`] ?? null],
      connectgaps: false,
    }));

    this.plotModalService.openPlot({
      data: traces,
      layout: {
        title: { text: this.trajectoryTitle },
        xaxis: { title: { text: 'Round' } },
        yaxis: { title: { text: this.trajectoryYAxis } },
        autosize: true,
      },
      config: { responsive: true, displaylogo: false },
    });
  }
}
