import { Component, inject, signal, output, OnDestroy } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { SplitPanel } from '../../common/split-panel/split-panel';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';
import { Table, TableConfig } from '../../common/table/table';

@Component({
  selector: 'app-count',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    SplitPanel,
    ...MATERIAL_IMPORTS],
  templateUrl: './count.html',
  styleUrl: './count.scss'
})
export class Count implements OnDestroy {

  tableConfig: TableConfig = {
    columns: [
      { key: 'id', label: 'id' },
      { key: 'rank', label: 'Rank' },
      { key: 'reads', label: 'Reads' },
      { key: 'rpm', label: 'RPU' },
      { key: 'length', label: 'Length' },
      { key: 'seqs', label: 'Sequence' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };
  
  private apiService = inject(ApiService);
  
  // Output events to emit to parent component
  resultsReady = output<any[]>();
  showReadsPerRankModal = output<{ data: any[], params: any }>();
  showSeqLengthModal = output<{ data: any, params: any }>();
  showAbundancePlotModal = output<{ data: any[], params: any }>();

  selectedFile: File | null = null;
  savedFileName: string = '';
  fileName: string = 'FASTQ or FASTA file';
  normalizeValue: string = '1e+06';
  returnReverseComplement: string = 'no';
  downloadFormat: string = 'fasta';
  uploadComplete: boolean = false;
  
  // Use signals for reactive state that affects the template
  isProcessing = signal(false);
  processedFileName = signal('');

  // Table data (temporary storage before emitting to parent)
  tableData: any[] = [];

  // Plot parameters
  minReadsToPlot: number = 10;
  maxRankToPlot: number = 100;
  
  // Reads per Rank plot customization
  adjustReadsPerRank: string = 'no';
  rprYMetric: string = 'reads';
  rprXAxis: string = 'Ranks of unique sequences';
  rprYAxis: string = 'Total reads per unique sequence';
  rprTitle: string = 'Read count for each rank';
  rprLineColor: string = '#87CEEB';
  
  // Sequence-length histogram customization
  adjustSeqLengthHistogram: string = 'no';
  histXAxis: string = 'Sequence length';
  histYAxis1: string = 'Unique sequences';
  histYAxis2: string = 'Read count';
  histTitle: string = 'Sequence-length histogram';
  histBarOutline: string = '#000000';
  histBarFill: string = '#87CEEB';
  histBarFill2: string = '#FFA500';
  
  // Abundance plot customization
  adjustAbundancePlot: string = 'no';
  useSingleton: string = 'yes';
  abundanceBreakpoints: string = '10,100,1000';
  abundanceXAxis: string = 'No. Reads';
  abundanceYAxis: string = 'Fraction of Population';
  abundancePlotTitle: string = 'Binned sequence abundance';
  abundanceBarOutline: string = '#000000';
  abundanceBarFill: string = '#87CEEB';
  abundanceColorLight: string = '#ADD8E6';
  abundanceColorDark: string = '#FF6B6B';
  
  // Statistics
  totalSequences = signal(0);
  uniqueSequences = signal(0);
  elapsedTime = signal('0.00');

  // Table preview is capped for very large result files (millions of rows would
  // freeze/crash the tab); totals above are still computed over the full file.
  readonly previewRowLimit = 50000;
  previewTruncated = signal(false);

  // Helper method for slider label formatting
  formatLabel(value: number): string {
    return `${value}`;
  }

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
    this.tableData = [];
    console.log('File selected:', result.fileName);
  }

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
  }

  onLoadResult(): void {
    if (!this.savedFileName) return;
    this.tableData = [];
    this.processedFileName.set(this.savedFileName);
    this.apiService.countPreview({ input_path: this.savedFileName, limit: this.previewRowLimit }).subscribe({
      next: (response) => this.applyPreviewResponse(response),
      error: (error) => alert(`Failed to load preview: ${error.error?.detail || error.message}`)
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
    this.tableData = [];

    const params = {
      input_path: this.savedFileName,
      reverseComplement: this.returnReverseComplement === 'yes',
      scaling_factor: parseFloat(this.normalizeValue),
      output_format: this.downloadFormat
    };

    console.log('Starting count with parameters:', params);

    // Chain the operations using RxJS operators
    this.apiService.count(params).pipe(
      tap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Count complete:', response.result);
        }
      }),
      switchMap(response => {
        // Automatically load a preview after successful count
        if (response.status === 'ok' && response.result) {
          return this.apiService.countPreview({ input_path: response.result, limit: this.previewRowLimit });
        }
        return of(null);
      }),
      tap(response => {
        if (response) this.applyPreviewResponse(response);
      }),
      catchError(error => {
        const errorMsg = error.error?.detail || error.message || 'Count failed';
        alert(`Count failed: ${errorMsg}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  // Preview rows + totals now come from the backend (/count-preview), computed
  // in a single streaming pass over the result file - the frontend never
  // downloads or parses the full file, so this scales to multi-GB results.
  private applyPreviewResponse(response: any): void {
    this.tableData = response.rows || [];
    this.uniqueSequences.set(response.total_records || 0);
    this.totalSequences.set(response.total_reads || 0);
    this.previewTruncated.set(!!response.truncated);
  }

  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) {
      console.warn('No file available for download. Please run count first.');
      return;
    }

    const link = document.createElement('a');
    link.href = this.apiService.getDownloadUrl(filename);
    link.download = filename;
    link.click();
  }

  // Plot methods - fetch aggregated data from the backend (computed over the
  // full result file, not just the capped table preview) and emit to parent.
  isLoadingReadsPerRank = signal(false);
  isLoadingSeqLengthHistogram = signal(false);
  isLoadingAbundancePlot = signal(false);

  openReadsPerRankPlot(): void {
    if (!this.processedFileName()) {
      alert('No data available for plotting. Please run the count process first.');
      return;
    }

    if (this.minReadsToPlot < 0) {
      alert('Minimum number of reads cannot be negative.');
      return;
    }

    if (this.maxRankToPlot < 1) {
      alert('Maximum rank must be at least 1.');
      return;
    }

    this.isLoadingReadsPerRank.set(true);
    this.apiService.countReadsPerRank({
      input_path: this.processedFileName(),
      min_reads: this.minReadsToPlot,
      max_rank: this.maxRankToPlot,
      metric: this.rprYMetric
    }).pipe(
      finalize(() => this.isLoadingReadsPerRank.set(false))
    ).subscribe({
      next: (response) => {
        const data = (response.data || []).map((d: any) => ({ rank: d.rank, reads: d.value }));
        if (data.length === 0) {
          alert('No data points match the specified criteria. Please adjust the min reads or max rank values.');
          return;
        }
        const defaultYLabel = this.rprYMetric === 'rpu' ? 'RPU per unique sequence' : 'Total reads per unique sequence';
        const params = {
          title: this.rprTitle,
          xAxisLabel: this.rprXAxis,
          yAxisLabel: this.adjustReadsPerRank === 'yes' ? this.rprYAxis : defaultYLabel,
          lineColor: this.rprLineColor
        };
        this.showReadsPerRankModal.emit({ data, params });
      },
      error: (error) => {
        alert(`Failed to load plot data: ${error.error?.detail || error.message}`);
      }
    });
  }

  openSeqLengthHistogram(): void {
    if (!this.processedFileName()) {
      alert('No data available for plotting. Please run the count process first.');
      return;
    }

    this.isLoadingSeqLengthHistogram.set(true);
    this.apiService.countSequenceLengthHistogram({
      input_path: this.processedFileName()
    }).pipe(
      finalize(() => this.isLoadingSeqLengthHistogram.set(false))
    ).subscribe({
      next: (response) => {
        const lengths: number[] = response.lengths || [];
        const data = {
          unique: lengths.map((l, i) => ({ length: l, count: response.unique[i] })),
          total: lengths.map((l, i) => ({ length: l, count: response.reads[i] }))
        };
        if (data.unique.length === 0 && data.total.length === 0) {
          alert('No sequence length data available for plotting.');
          return;
        }
        const outliers = response.excluded_outliers || [];
        if (outliers.length > 0) {
          const preview = outliers.slice(0, 5).map((o: any) => `${o.length} nt (${o.unique} seq)`).join(', ');
          alert(
            `${outliers.length} sequence length(s) far outside the normal range were excluded from this chart ` +
            `so the real distribution stays readable: ${preview}${outliers.length > 5 ? ', ...' : ''}. ` +
            `These are likely malformed records rather than real sequences - the full result file (via Download) still contains them.`
          );
        }
        const params = {
          title: this.histTitle,
          xAxisLabel: this.histXAxis,
          yAxis1Label: this.histYAxis1,
          yAxis2Label: this.histYAxis2,
          barOutline: this.histBarOutline,
          barFill: this.histBarFill,
          barFill2: this.histBarFill2
        };
        this.showSeqLengthModal.emit({ data, params });
      },
      error: (error) => {
        alert(`Failed to load plot data: ${error.error?.detail || error.message}`);
      }
    });
  }

  openAbundancePlot(): void {
    if (!this.processedFileName()) {
      alert('No data available for plotting. Please run the count process first.');
      return;
    }

    const breaks = this.abundanceBreakpoints.split(',').map(b => parseInt(b.trim())).filter(b => !isNaN(b));

    if (breaks.length === 0) {
      alert('Invalid abundance breakpoints. Please enter comma-separated numbers (e.g., 10,100,1000).');
      return;
    }

    if (breaks.some(b => b <= 0)) {
      alert('Abundance breakpoints must be positive numbers.');
      return;
    }

    this.isLoadingAbundancePlot.set(true);
    this.apiService.countAbundance({
      input_path: this.processedFileName(),
      breakpoints: breaks,
      use_singleton: this.useSingleton === 'yes'
    }).pipe(
      finalize(() => this.isLoadingAbundancePlot.set(false))
    ).subscribe({
      next: (response) => {
        const data = (response.data || []).map((d: any) => ({
          bin: d.bin, fraction: d.fraction, uniqueCount: d.unique_count
        }));
        if (data.length === 0) {
          alert('No abundance data available for plotting with the specified breakpoints.');
          return;
        }
        const params = {
          title: this.abundancePlotTitle,
          xAxisLabel: this.abundanceXAxis,
          yAxisLabel: this.abundanceYAxis,
          barOutline: this.abundanceBarOutline,
          barFill: this.abundanceBarFill,
          colorLight: this.abundanceColorLight,
          colorDark: this.abundanceColorDark
        };
        this.showAbundancePlotModal.emit({ data, params });
      },
      error: (error) => {
        alert(`Failed to load plot data: ${error.error?.detail || error.message}`);
      }
    });
  }

}
