import { Component, inject, signal, ChangeDetectorRef, OnDestroy } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { ApiService } from '../../../shared/api.service';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { Table, TableConfig } from '../../common/table/table';
import { SplitPanel } from '../../common/split-panel/split-panel';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';

interface UploadedFile {
  savedFileName: string;
  displayName: string;
}

@Component({
  selector: 'app-diff-analysis',
  imports: [CommonModule, FormsModule, Upload, Table, SplitPanel, ...MATERIAL_IMPORTS],
  templateUrl: './diff-analysis.html',
  styleUrl: './diff-analysis.scss',
  standalone: true
})
export class DiffAnalysis implements OnDestroy {
  private apiService = inject(ApiService);
  private plotModalService = inject(PlotModalService);
  private cdr = inject(ChangeDetectorRef);

  isProcessing = signal(false);
  processedFileName = signal('');

  filesCond1: UploadedFile[] = [];
  filesCond2: UploadedFile[] = [];

  lfcCutoff: number = 1.0;
  downloadFormat: string = 'csv';
  diffAnalysisData: any[] = [];

  tableConfig: TableConfig = {
    columns: [
      { key: 'Sequence', label: 'Sequence' },
      { key: 'logFC', label: 'logFC', exact_match: true },
      { key: 'logCPM', label: 'logCPM', exact_match: true },
      { key: 'PValue', label: 'PValue', exact_match: true },
      { key: 'PClass', label: 'PClass' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  adjustEdgeRPlot = 'no';
  edgeRPlotXAxis = 'logCPM';
  edgeRPlotYAxis = 'logFC';
  edgeRPlotTitle = 'edgeR results';
  edgeRPlotSigColor = '#FF0000';
  edgeRPlotInsigColor = '#000000';

  onAdjustEdgeRPlotChange(): void {
    if (this.adjustEdgeRPlot === 'no') {
      this.edgeRPlotXAxis = 'logCPM';
      this.edgeRPlotYAxis = 'logFC';
      this.edgeRPlotTitle = 'edgeR results';
      this.edgeRPlotSigColor = '#FF0000';
      this.edgeRPlotInsigColor = '#000000';
    }
  }

  // ========================================================================
  // UPLOAD HANDLERS
  // ========================================================================

  onCond1Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.filesCond1.push({ savedFileName: result.savedFileName, displayName: result.fileName });
    }
  }

  onCond2Complete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.filesCond2.push({ savedFileName: result.savedFileName, displayName: result.fileName });
    }
  }

  removeCond1File(file: UploadedFile): void {
    this.apiService.deleteFile(file.savedFileName).subscribe();
    this.filesCond1 = this.filesCond1.filter(f => f !== file);
  }

  removeCond2File(file: UploadedFile): void {
    this.apiService.deleteFile(file.savedFileName).subscribe();
    this.filesCond2 = this.filesCond2.filter(f => f !== file);
  }

  get canRun(): boolean {
    return this.filesCond1.length >=1 && this.filesCond2.length >= 1 && !this.isProcessing();
  }

  // ========================================================================
  // ANALYSIS
  // ========================================================================

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
  }

  ngOnDestroy(): void {
    [...this.filesCond1, ...this.filesCond2].forEach(f => this.apiService.deleteFile(f.savedFileName).subscribe());
    if (this.processedFileName()) this.apiService.deleteFile(this.processedFileName()).subscribe();
  }

  onStart(): void {
    this.isProcessing.set(true);
    this.diffAnalysisData = [];

    this.apiService.differentialAnalysis({
      cond1_paths: this.filesCond1.map(f => f.savedFileName),
      cond2_paths: this.filesCond2.map(f => f.savedFileName),
      lfc_cutoff: this.lfcCutoff,
      output_format: this.downloadFormat
    }).pipe(
      tap(response => this.processedFileName.set(response.result)),
      switchMap(response => {
        if (!response.result) throw new Error('No result file returned');
        return this.apiService.downloadFile(response.result);
      }),
      tap(blob => this.parseFileBlob(blob)),
      catchError(error => {
        alert(`Differential analysis failed: ${error.error?.detail || error.message}`);
        return of(null);
      }),
      finalize(() => { this.isProcessing.set(false); this.cdr.detectChanges(); })
    ).subscribe();
  }

  private parseFileBlob(blob: Blob): void {
    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      this.diffAnalysisData = [];
      const lines = content.split('\n').filter(l => l.trim());
      if (lines.length === 0) return;
      const headers = lines[0].split(',').map(h => h.trim());
      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',').map(v => v.trim());
        const row: any = {};
        headers.forEach((header, idx) => {
          const value = values[idx];
          if (header === 'PValue') {
            row[header] = (!value || value === '' || value.toLowerCase() === 'nan') ? 'N/A' : parseFloat(value);
          } else if (['logFC', 'logCPM'].includes(header)) {
            row[header] = value ? parseFloat(value) : 0;
          } else {
            row[header] = value || '';
          }
        });
        this.diffAnalysisData.push(row);
      }
      this.cdr.detectChanges();
    };
    reader.readAsText(blob);
  }

  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) return;
    this.apiService.downloadFile(filename).subscribe({
      next: (blob: Blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        window.URL.revokeObjectURL(url);
      }
    });
  }

  // ========================================================================
  // EDGER PLOT
  // ========================================================================

  edgeRPlot(): void {
    if (this.diffAnalysisData.length === 0) return;

    const logCPM = this.diffAnalysisData.map(r => r['logCPM']);
    const logFC = this.diffAnalysisData.map(r => r['logFC']);
    const pClass = this.diffAnalysisData.map(r => r['PClass']);
    const sequences = this.diffAnalysisData.map(r => r['Sequence'] || '');

    const sigIdx = pClass.map((p, i) => p === 'Sig.' ? i : -1).filter(i => i >= 0);
    const insigIdx = pClass.map((p, i) => p === 'Not Sig.' ? i : -1).filter(i => i >= 0);
    const traces: any[] = [];

    if (insigIdx.length > 0) {
      traces.push({
        x: insigIdx.map(i => logCPM[i]), y: insigIdx.map(i => logFC[i]),
        mode: 'markers', type: 'scatter',
        marker: { color: this.edgeRPlotInsigColor, size: 6, opacity: 0.5 },
        text: insigIdx.map(i => sequences[i]),
        hovertemplate: '%{text}<br>logCPM: %{x:.3f}<br>logFC: %{y:.3f}<extra></extra>',
        name: 'Not Significant'
      });
    }
    if (sigIdx.length > 0) {
      traces.push({
        x: sigIdx.map(i => logCPM[i]), y: sigIdx.map(i => logFC[i]),
        mode: 'markers', type: 'scatter',
        marker: { color: this.edgeRPlotSigColor, size: 6, opacity: 0.7 },
        text: sigIdx.map(i => sequences[i]),
        hovertemplate: '%{text}<br>logCPM: %{x:.3f}<br>logFC: %{y:.3f}<extra></extra>',
        name: 'Significant'
      });
    }

    this.plotModalService.openPlot({
      data: traces,
      layout: {
        title: { text: this.edgeRPlotTitle },
        xaxis: { title: { text: this.edgeRPlotXAxis } },
        yaxis: { title: { text: this.edgeRPlotYAxis } },
        autosize: true, hovermode: 'closest', showlegend: true
      },
      config: { responsive: true }
    });
  }
}
