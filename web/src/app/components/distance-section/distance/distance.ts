import { Component, inject, signal, ChangeDetectorRef, OnDestroy } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { FileService, ColumnName } from '../../../shared/file-service';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { Table, TableConfig } from '../../common/table/table';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';

@Component({
  selector: 'app-distance',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './distance.html',
  styleUrl: './distance.scss',
  standalone: true
})
export class Distance implements OnDestroy {
  // Services
  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private plotModalService = inject(PlotModalService);
  private cdr = inject(ChangeDetectorRef);

  // Signals
  isProcessing = signal(false);
  processedFileName = signal('');

  // File handling
  savedFileName: string = '';
  uploadComplete: boolean = false;

  // Query sequence
  querySequence: string = '';

  // Download format
  downloadFormat: string = 'csv';

  // Data State
  distanceData: any[] = [];

  // Table configuration
  tableConfig: TableConfig = {
    columns: [
      { key: ColumnName.ID, label: 'ID' },
      { key: ColumnName.RANK, label: 'Rank', exact_match: true },
      { key: ColumnName.READS, label: 'Reads', exact_match: true },
      { key: ColumnName.RPU, label: 'RPU', exact_match: true },
      { key: ColumnName.SEQUENCES, label: 'Sequences' },
      { key: ColumnName.DISTANCE, label: 'Distance', exact_match: true }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  // Distance histogram customization
  adjustHistogram: string = 'no';
  histXAxis: string = 'Distance from query';
  histYAxis1: string = 'Unique sequences';
  histYAxis2: string = 'Read count';
  histTitle: string = 'Distance histograms';
  histBarOutline: string = '#000000';
  histBarFill: string = '#87CEEB';

  // Watch for changes in adjustment toggle
  onAdjustHistogramChange(): void {
    if (this.adjustHistogram === 'no') {
      this.histXAxis = 'Distance from query';
      this.histYAxis1 = 'Unique sequences';
      this.histYAxis2 = 'Read count';
      this.histTitle = 'Distance histograms';
      this.histBarOutline = '#000000';
      this.histBarFill = '#87CEEB';
    }
  }

  // File handlers
  onFileSelected(result: FileUploadResult): void {
    if (this.savedFileName) {
      this.apiService.deleteFile(this.savedFileName).subscribe();
    }
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
    this.savedFileName = '';
    this.uploadComplete = false;
    this.processedFileName.set('');
    this.distanceData = [];
  }

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
  }

  onLoadResult(): void {
    if (!this.savedFileName) return;
    this.distanceData = [];
    this.apiService.downloadFile(this.savedFileName).pipe(
      switchMap(blob => this.fileService.parseClusterFile(blob, this.savedFileName)),
      tap(parsedData => {
        this.distanceData = parsedData;
        this.processedFileName.set(this.savedFileName);
      })
    ).subscribe();
  }

  ngOnDestroy(): void {
    if (this.savedFileName) {
      this.apiService.deleteFile(this.savedFileName).subscribe();
    }
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
  }

  onFileUploadComplete(result: FileUploadResult): void {
    this.savedFileName = result.savedFileName || '';
    this.uploadComplete = true;
  }

  onStart(): void {
    if (!this.uploadComplete) {
      alert('Please upload a FASTA file before starting.');
      return;
    }

    if (!this.savedFileName) {
      alert('File upload incomplete. Please try again.');
      return;
    }

    if (!this.querySequence || this.querySequence.trim() === '') {
      alert('Please enter a query sequence.');
      return;
    }

    // Validate query sequence (alphanumeric only)
    const alphanumericRegex = /^[a-zA-Z0-9]+$/;
    if (!alphanumericRegex.test(this.querySequence.trim())) {
      alert('Query sequence must be alphanumeric (letters and numbers only).');
      return;
    }

    this.isProcessing.set(true);
    this.distanceData = [];

    const params = {
      input_path: this.savedFileName,
      query_sequence: this.querySequence.trim(),
      output_format: this.downloadFormat
    };

    this.apiService.sequenceDistance(params).pipe(
      switchMap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Distance completed:', response.result);
          
          return this.apiService.downloadFile(response.result).pipe(
            switchMap(blob => 
              this.fileService.parseClusterFile(blob, response.result)
            ),
            tap(parsedData => {
              this.distanceData = parsedData;
              this.cdr.detectChanges();
            })
          );
        }
        return of(null);
      }),
      catchError(error => {
        const errorMessage = error.error?.detail || error.message || 'Unknown error';
        alert(`Distance calculation failed: ${errorMessage}`);
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
      console.warn('No file available for download. Please run distance calculation first.');
      return;
    }

    console.log('Downloading file:', filename);
    this.fileService.downloadFile(filename);
  }

  // ========================================================================
  // PLOTTING
  // ========================================================================

  async distanceHistogram(): Promise<void> {
    if (this.distanceData.length === 0) {
      alert('No distance data available for plotting. Please run the analysis first.');
      return;
    }

    // Extract distance values
    const distances = this.distanceData
      .map(row => row[ColumnName.DISTANCE])
      .filter(val => val !== null && val !== undefined && isFinite(val));

    if (distances.length === 0) {
      alert('No valid distance values found for plotting.');
      return;
    }

    // Calculate histogram bins for unique sequences
    const maxDistance = Math.max(...distances);
    const binSize = 1;
    const bins: { [key: number]: number } = {};
    
    distances.forEach(d => {
      bins[d] = (bins[d] || 0) + 1;
    });

    const uniqueX: number[] = [];
    const uniqueY: number[] = [];
    for (let i = 0; i <= maxDistance; i++) {
      uniqueX.push(i);
      uniqueY.push(bins[i] || 0);
    }

    // Calculate histogram for read counts
    const readBins: { [key: number]: number } = {};
    this.distanceData.forEach(row => {
      const dist = row[ColumnName.DISTANCE];
      const reads = row[ColumnName.READS] || 0;
      readBins[dist] = (readBins[dist] || 0) + reads;
    });

    const readsX: number[] = [];
    const readsY: number[] = [];
    for (let i = 0; i <= maxDistance; i++) {
      readsX.push(i);
      readsY.push(readBins[i] || 0);
    }

    // Create first trace for unique sequences
    const trace1 = {
      x: uniqueX,
      y: uniqueY,
      type: 'bar',
      marker: {
        color: this.histBarFill,
        line: {
          color: this.histBarOutline,
          width: 1
        }
      },
      name: 'Unique sequences'
    };

    // Create second trace for read counts
    const trace2 = {
      x: readsX,
      y: readsY,
      type: 'bar',
      marker: {
        color: this.histBarFill,
        line: {
          color: this.histBarOutline,
          width: 1
        }
      },
      name: 'Read count'
    };

    // Create subplots layout
    const layout = {
      title: {
        text: this.histTitle,
        font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
      },
      grid: { rows: 2, columns: 1, pattern: 'independent' },
      xaxis: {
        title: {
          text: '', // No title on top plot to avoid confusion
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0',
        anchor: 'y'
      },
      yaxis: {
        title: {
          text: this.histYAxis1,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0',
        domain: [0.55, 1]
      },
      xaxis2: {
        title: {
          text: this.histXAxis,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0',
        anchor: 'y2'
      },
      yaxis2: {
        title: {
          text: this.histYAxis2,
          font: { size: 14, family: 'Arial, sans-serif', weight: 'bold' }
        },
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        showgrid: true,
        gridcolor: '#e0e0e0',
        domain: [0, 0.45],
        anchor: 'x2'
      },
      autosize: true,
      height: 700,
      margin: { t: 60, b: 70, l: 90, r: 50 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white',
      showlegend: false
    };

    // Assign traces to different y-axes
    const trace1WithYaxis = { ...trace1, xaxis: 'x', yaxis: 'y' };
    const trace2WithYaxis = { ...trace2, xaxis: 'x2', yaxis: 'y2' };

    this.plotModalService.openPlot({
      data: [trace1WithYaxis, trace2WithYaxis],
      layout: layout,
      config: { responsive: true }
    });
  }
}
