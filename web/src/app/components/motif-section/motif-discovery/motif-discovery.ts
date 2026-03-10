import { Component, inject, signal, ChangeDetectorRef, NgZone, OnDestroy } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';
import { Table, TableConfig } from '../../common/table/table';

@Component({
  selector: 'app-motif-discovery',
  imports: [ 
    CommonModule,
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './motif-discovery.html',
  styleUrl: './motif-discovery.scss'
})
export class MotifDiscovery implements OnDestroy {

  tableConfig: TableConfig = {
    columns: [
      { key: 'motif', label: 'Motif' },
      { key: 'p', label: 'P-value' },
      { key: 'zz', label: 'ZZ Score' },
      { key: 'motif_length', label: 'Length' },
      { key: 'rank', label: 'Rank' },
      { key: 'seq_count', label: 'Seq Count' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  private apiService = inject(ApiService);
  private cdr = inject(ChangeDetectorRef);
  private plotModalService = inject(PlotModalService);
  private ngZone = inject(NgZone);

  selectedFile: File | null = null;
  savedFileName: string = '';
  fileName: string = 'FASTA file';
  minReads: number = 10;
  minLength: number = 5;
  maxLength: number = 10;
  uploadComplete: boolean = false;
  alphabet: string = 'dna';
  
  // Plot customization options
  showPlotCustomization: string = 'no';
  plotXaxis: string = 'Rank by normalized z-score';
  plotYaxis: string = '-log10(p)';
  plotLegend: string = 'Length';
  plotTitle: string = 'Over-enriched strings';
  plotPalette: string = 'magma';
  
  // Default plot values
  private readonly defaultPlotValues = {
    plotXaxis: 'Rank by normalized z-score',
    plotYaxis: '-log10(p)',
    plotLegend: 'Length',
    plotTitle: 'Over-enriched strings',
    plotPalette: 'magma'
  };
  
  paletteOptions = ['magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo'];

  // Use signals for reactive state
  isProcessing = signal(false);
  processedFileName = signal('');

  // Table data
  tableData: any[] = [];

  onPlotCustomizationChange(): void {
    if (this.showPlotCustomization === 'no') {
      // Reset to defaults when user selects "No"
      this.plotXaxis = this.defaultPlotValues.plotXaxis;
      this.plotYaxis = this.defaultPlotValues.plotYaxis;
      this.plotLegend = this.defaultPlotValues.plotLegend;
      this.plotTitle = this.defaultPlotValues.plotTitle;
      this.plotPalette = this.defaultPlotValues.plotPalette;
    }
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
      alert('Please upload a file first!');
      return;
    }

    if (this.minLength > this.maxLength) {
      alert('Minimum length must be less than maximum length!');
      return;
    }

    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.tableData = [];

    const params = {
      input_path: this.savedFileName,
      min_reads: this.minReads,
      length_range: [this.minLength, this.maxLength],
      output_format: 'csv',
      alphabet: this.alphabet
    };

    console.log('Starting motif discovery with parameters:', params);

    // Chain the operations using RxJS operators
    this.apiService.motifDiscovery(params).pipe(
      tap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Motif discovery complete:', response.result);
        }
      }),
      switchMap(response => {
        // Automatically load results after successful discovery
        if (response.status === 'ok' && response.result) {
          return this.apiService.downloadFile(response.result).pipe(
            tap(blob => this.parseFileBlob(blob, response.result))
          );
        }
        return of(null);
      }),
      catchError(error => {
        const errorMsg = error.error?.detail || 'Motif discovery failed';
        console.error('Motif discovery error:', errorMsg);
        alert(`Error: ${errorMsg}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  parseFileBlob(blob: Blob, filename: string): void {
    const reader = new FileReader();
    reader.onload = (e: any) => {
      const text = e.target.result;
      this.ngZone.run(() => {
        this.parseResultFile(text, filename);
        this.cdr.detectChanges();
      });
    };
    reader.readAsText(blob);
  }

  parseResultFile(content: string, filename: string): void {
    const isCsv = filename.endsWith('.csv');
    
    this.tableData = [];
    
    if (isCsv) {
      // Parse CSV
      const lines = content.split('\n').filter(line => line.trim());
      if (lines.length === 0) return;
      
      const headers = lines[0].split(',').map(h => h.trim());
      
      // Find column indices (case-insensitive)
      const motifIdx = headers.findIndex(h => h.toLowerCase() === 'motif');
      const pIdx = headers.findIndex(h => h.toLowerCase() === 'p');
      const zzIdx = headers.findIndex(h => h.toLowerCase() === 'zz');
      const lengthIdx = headers.findIndex(h => h.toLowerCase() === 'motif_length');
      const rankIdx = headers.findIndex(h => h.toLowerCase() === 'rank');
      const seqCountIdx = headers.findIndex(h => h.toLowerCase() === 'seqcount');
      
      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',');
        if (values.length > 0 && motifIdx >= 0) {
          this.tableData.push({
            motif: values[motifIdx] || '',
            p: pIdx >= 0 ? parseFloat(values[pIdx]) : 0,
            zz: zzIdx >= 0 ? parseFloat(values[zzIdx]) : 0,
            motif_length: lengthIdx >= 0 ? parseInt(values[lengthIdx]) : 0,
            rank: rankIdx >= 0 ? parseInt(values[rankIdx]) : 0,
            seq_count: seqCountIdx >= 0 ? parseInt(values[seqCountIdx]) : 0
          });
        }
      }
    }
    
    console.log('Parsed motif discovery results:', this.tableData.length, 'motifs');
  }

  onDownload(): void {
    if (!this.processedFileName()) {
      console.warn('No processed file available for download');
      return;
    }

    this.apiService.downloadFile(this.processedFileName()).subscribe({
      next: (blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = this.processedFileName();
        link.click();
        window.URL.revokeObjectURL(url);
        console.log('File downloaded:', this.processedFileName());
      },
      error: (error) => {
        console.error('Download error:', error);
        alert('Failed to download file');
      }
    });
  }

  onShowPlot(): void {
    if (this.tableData.length === 0) {
      alert('No data available to plot. Please run motif discovery first.');
      return;
    }

    // Get the color palette name
    const colorPalette = this.plotPalette;
    
    // Prepare data for Plotly bubble plot
    const plotData = [{
      x: this.tableData.map(d => d.rank),
      y: this.tableData.map(d => -Math.log10(d.p)),
      mode: 'markers',
      marker: {
        size: this.tableData.map(d => d.motif_length * 2), // Scale size for visibility
        color: this.tableData.map(d => d.motif_length),
        colorscale: this.getColorscale(colorPalette),
        colorbar: {
          title: { text: this.plotLegend }
        },
        line: {
          color: 'black',
          width: 0.5
        }
      },
      text: this.tableData.map(d => 
        `Motif: ${d.motif}<br>Rank: ${d.rank}<br>-log10(p-value): ${(-Math.log10(d.p)).toFixed(2)}<br>Motif length: ${d.motif_length}`
      ),
      hoverinfo: 'text',
      type: 'scatter'
    }];

    const layout = {
      title: {
        text: this.plotTitle,
        font: { size: 16 }
      },
      xaxis: {
        title: {
          text: this.plotXaxis,
          font: { size: 14 }
        }
      },
      yaxis: {
        title: {
          text: this.plotYaxis,
          font: { size: 14 }
        }
      },
      hovermode: 'closest',
      showlegend: false
    };

    this.plotModalService.openPlot({
      data: plotData,
      layout: layout
    });
  }

  // Map palette names to Plotly colorscales
  private getColorscale(paletteName: string): string {
    const colorscaleMap: { [key: string]: string } = {
      'magma': 'Hot',
      'inferno': 'Hot',
      'plasma': 'Portland',
      'viridis': 'Viridis',
      'cividis': 'Cividis',
      'rocket': 'Jet',
      'mako': 'Blues',
      'turbo': 'Rainbow'
    };
    
    return colorscaleMap[paletteName] || 'Viridis';
  }
}
