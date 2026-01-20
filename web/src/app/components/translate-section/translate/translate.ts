import { Component, inject, signal } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { ColumnName, FileService } from '../../../shared/file-service';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';
import { Table, TableConfig } from '../../common/table/table';
import { PlotModalService } from '../../../shared/plot-modal.service';

export interface CodonChange {
  Codon: string;
  Translation: string;
}

@Component({
  selector: 'app-translate',
  standalone: true,
  imports: [ 
    CommonModule, 
    FormsModule,
    Upload,
    Table,
    ...MATERIAL_IMPORTS],
  templateUrl: './translate.html',
  styleUrl: './translate.scss'
})
export class Translate {

  // Services
  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private plotModalService = inject(PlotModalService);

  // Signals
  isProcessing = signal(false);
  processedFileName = signal('');

  // File handling
  selectedFile: File | null = null;
  savedFileName: string = '';
  uploadComplete: boolean = false;

  // Data State
  translateData: any[] = [];

  tableConfig: TableConfig = {
    columns: [
      { key: ColumnName.ID, label: 'ID' },
      { key: ColumnName.RANK, label: 'Rank', exact_match: true },
      { key: ColumnName.READS, label: 'Reads', exact_match: true },
      { key: ColumnName.RPU, label: 'RPU', exact_match: true },
      { key: ColumnName.UNIQUE_NTS, label: 'Unique Nts', exact_match: true },
      { key: ColumnName.LENGTH, label: 'Length', exact_match: true },
      { key: ColumnName.SEQUENCES, label: 'Sequence' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };
  
  // Translation parameters
  orf: number = 1;
  converge: string = 'yes';
  translateSelection: string = 'Standard';
  downloadFormat: string = 'fasta';
  
  // Custom translation table
  nonstandardTranslations: string = 'no';
  translateInputChangesCodons: string = '';
  translateInputChangesOutputs: string = '';

  // Plot parameters
  minReadsToPlot: number = 10;
  maxRankToPlot: number = 100;
  
  // Reads per Rank plot customization
  adjustReadsPerRank: string = 'no';
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

  // Available genetic codes
  geneticCodes: string[] = [
    'Standard',
    'Vertebrate mitochondrial',
    'Yeast mitochondrial',
    'Mold, protozoan, and coelenterate mitochondrial + Mycoplasma / Spiroplasma',
    'Invertebrate mitochondrial',
    'Ciliate, dasycladacean and Hexamita nuclear',
    'Echinoderm and flatworm mitochondrial',
    'Euplotid nuclear',
    'Alternative yeast nuclear',
    'Ascidian mitochondrial',
    'Alternative flatworm mitochondrial',
    'Blepharisma nuclear',
    'Chlorophycean mitochondrial',
    'Trematode mitochondrial',
    'Scenedesmus obliquus mitochondrial',
    'Pterobranchia mitochondrial'
  ];

  // ========================================================================
  // FILE HANDLING
  // ========================================================================

  // Reset plot customizations to defaults
  resetReadsPerRankDefaults(): void {
    this.rprXAxis = 'Ranks of unique sequences';
    this.rprYAxis = 'Total reads per unique sequence';
    this.rprTitle = 'Read count for each rank';
    this.rprLineColor = '#87CEEB';
  }

  resetSeqLengthHistogramDefaults(): void {
    this.histXAxis = 'Sequence length';
    this.histYAxis1 = 'Unique sequences';
    this.histYAxis2 = 'Read count';
    this.histTitle = 'Sequence-length histogram';
    this.histBarOutline = '#000000';
    this.histBarFill = '#87CEEB';
    this.histBarFill2 = '#FFA500';
  }

  // Watch for changes in adjustment toggles
  onAdjustReadsPerRankChange(): void {
    if (this.adjustReadsPerRank === 'no') {
      this.resetReadsPerRankDefaults();
    }
  }

  onAdjustSeqLengthHistogramChange(): void {
    if (this.adjustSeqLengthHistogram === 'no') {
      this.resetSeqLengthHistogramDefaults();
    }
  }

  onFileSelected(result: FileUploadResult): void {
    this.selectedFile = result.file;
    this.processedFileName.set('');
    this.translateData = [];
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

  onStart() {
    if (!this.uploadComplete || !this.savedFileName) {
      console.error('No file uploaded');
      return;
    }

    this.isProcessing.set(true);
    this.translateData = [];

    // Parse custom codon changes if provided
    let inputChanges: CodonChange[] | null = null;
    if (this.nonstandardTranslations === 'yes' && 
        this.translateInputChangesCodons && 
        this.translateInputChangesOutputs) {
      
      const codons = this.translateInputChangesCodons
        .split(',')
        .map(c => c.trim())
        .filter(c => c.length > 0);
      
      const translations = this.translateInputChangesOutputs
        .split(',')
        .map(t => t.trim())
        .filter(t => t.length > 0);
      
      if (codons.length === translations.length && codons.length > 0) {
        inputChanges = codons.map((codon, i) => ({
          Codon: codon,
          Translation: translations[i]
        }));
      }
    }

    const params = {
      input_path: this.savedFileName,
      orf: this.orf,
      converge: this.converge === 'yes',
      input_changes: inputChanges,
      translate_selection: this.translateSelection,
      output_format: this.downloadFormat
    };

    this.apiService.translate(params).pipe(
      switchMap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Translation completed:', response.result);
          
          // Chain download and parsing
          return this.apiService.downloadFile(response.result).pipe(
            switchMap(blob => this.fileService.parseClusterFile(blob, response.result)),
            tap(parsedData => {
              this.translateData = parsedData;
            })
          );
        }
        return of(null);
      }),
      catchError(error => {
        console.error('Translation error:', error);
        alert('Translation failed!');
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  // ========================================================================
  // PLOTTING
  // ========================================================================

  async readsPerRankPlot() {
    if (this.translateData.length === 0) {
      console.error('No data available for plotting');
      return;
    }

    // Filter data based on min reads and max rank
    const filteredData = this.translateData
      .filter((row: any) => row[ColumnName.READS] >= this.minReadsToPlot)
      .slice(0, this.maxRankToPlot);

    if (filteredData.length === 0) {
      alert(`No data to plot! All sequences have reads < ${this.minReadsToPlot}. Try lowering the minimum reads threshold.`);
      return;
    }

    // Prepare plot data
    const xData = filteredData.map((row: any) => row[ColumnName.RANK]);
    const yData = filteredData.map((row: any) => row[ColumnName.READS]);

    const trace = {
      x: xData,
      y: yData,
      type: 'scatter',
      mode: 'lines+markers',
      line: { color: this.rprLineColor, width: 2 },
      marker: { size: 6, color: this.rprLineColor }
    };

    const layout = {
      title: {
        text: this.rprTitle,
        font: { size: 18 }
      },
      xaxis: {
        title: {
          text: this.rprXAxis,
          font: { size: 14 }
        },
        showgrid: true,
        gridcolor: '#e0e0e0',
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        mirror: true
      },
      yaxis: {
        title: {
          text: this.rprYAxis,
          font: { size: 14 }
        },
        showgrid: true,
        gridcolor: '#e0e0e0',
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        mirror: true
      },
      autosize: true,
      height: 500,
      margin: { t: 60, b: 70, l: 90, r: 40 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white'
    };

    // Use shared plot modal service
    this.plotModalService.openPlot({
      data: [trace],
      layout: layout,
      config: { responsive: true }
    });
  }

  async seqLengthHistogramPlot() {
    if (this.translateData.length === 0) {
      console.error('No data available for plotting');
      return;
    }

    // Group data by sequence length
    const lengthCounts = new Map<number, { unique: number, total: number }>();
    
    this.translateData.forEach((row: any) => {
      const length = row[ColumnName.LENGTH] || row[ColumnName.SEQUENCES]?.length || 0;
      const reads = row[ColumnName.READS] || 1;
      
      if (!lengthCounts.has(length)) {
        lengthCounts.set(length, { unique: 0, total: 0 });
      }
      
      const counts = lengthCounts.get(length)!;
      counts.unique += 1;
      counts.total += reads;
    });

    // Convert to sorted arrays
    const sortedEntries = Array.from(lengthCounts.entries()).sort((a, b) => a[0] - b[0]);
    const lengths = sortedEntries.map(([length]) => length);
    const uniqueCounts = sortedEntries.map(([, counts]) => counts.unique);
    const totalCounts = sortedEntries.map(([, counts]) => counts.total);

    const trace1 = {
      x: lengths,
      y: uniqueCounts,
      name: this.histYAxis1,
      type: 'bar',
      marker: { 
        color: this.histBarFill,
        line: { color: this.histBarOutline, width: 1 }
      },
      yaxis: 'y'
    };

    const trace2 = {
      x: lengths,
      y: totalCounts,
      name: this.histYAxis2,
      type: 'bar',
      marker: { 
        color: this.histBarFill2,
        line: { color: this.histBarOutline, width: 1 }
      },
      yaxis: 'y2'
    };

    const layout = {
      title: {
        text: this.histTitle,
        font: { size: 18 }
      },
      xaxis: {
        title: {
          text: this.histXAxis,
          font: { size: 14 }
        },
        showgrid: true,
        gridcolor: '#e0e0e0',
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        mirror: true
      },
      yaxis: {
        title: {
          text: this.histYAxis1,
          font: { size: 14 }
        },
        side: 'left',
        showgrid: true,
        gridcolor: '#e0e0e0',
        showline: true,
        linewidth: 2,
        linecolor: 'black'
      },
      yaxis2: {
        title: {
          text: this.histYAxis2,
          font: { size: 14 }
        },
        side: 'right',
        overlaying: 'y',
        showgrid: false,
        showline: true,
        linewidth: 2,
        linecolor: 'black'
      },
      autosize: true,
      height: 500,
      barmode: 'group',
      margin: { t: 60, b: 70, l: 90, r: 90 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white'
    };

    // Use shared plot modal service
    this.plotModalService.openPlot({
      data: [trace1, trace2],
      layout: layout,
      config: { responsive: true }
    });
  }
}
