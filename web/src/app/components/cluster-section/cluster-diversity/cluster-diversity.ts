import { Component, inject, signal, output } from '@angular/core';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { ApiService } from '../../../shared/api.service';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';
import { KmerAnalysisRequest } from '../../../shared/kmer-analysis.types';

export interface DiversityResultsEvent {
  data: any[];
  inputFile?: string;
}

@Component({
  selector: 'app-cluster-diversity',
  imports: [   
    CommonModule,
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS],
  templateUrl: './cluster-diversity.html',
  styleUrl: './cluster-diversity.scss',
})
export class ClusterDiversity {
  private apiService = inject(ApiService);

  // Output events for plot modals
  showMetaplotModal = output<{ data: any, params: any }>();
  showKmerPlotModal = output<{ data: any, params: any }>();
  
  // Output event for results ready
  resultsReady = output<DiversityResultsEvent>();

  // UI state
  adjustClusterMetadata: string = 'no';
  clusterMetaXAxis: string = "Cluster";
  clusterMetaYAxis1: string = "Seq. count";
  clusterMetaColor1: string = "#0000ff";
  clusterMetaYAxis2: string = "Read count";
  clusterMetaColor2: string = "#ffa500";
  clusterMetaYAxis3: string = "Avg. LED";
  clusterMetaColor3: string = "#228b22";
  clusterMetaPlotTitle: string = "Cluster metaplots";

  kmerSize: string = "3";
  plotType: string = "UMAP";

  adjustKmerPlot: string = 'no';
  kmerPlotXAxis: string = "Dim1";
  kmerPlotYAxis: string = "Dim2";
  kmerPlotLegendTitle: string = "Cluster";
  kmerPlotTitle: string = "Cluster k-mer plot";
  kmerPlotColorPalette: string = "Dark2";

  // File upload state
  selectedFile: File | null = null;
  savedFileName: string = '';
  uploadComplete: boolean = false;
  
  // Processing state
  isProcessing = signal(false);
  processedFileName = signal('');
  outputFormat: string = 'csv';

  // Cluster selection for k-mer plot
  availableClusters: number[] = [];
  selectedClusters: number[] = [];

  // Table data storage
  diversityData: any[] = [];
  displayedColumns: string[] = ['cluster', 'totalSequences', 'totalReads', 'totalRPU', 'averageLED', 'SID', 'seedID', 'seedSequence'];

  onFileSelected(result: FileUploadResult): void {
    this.selectedFile = result.file;
    this.processedFileName.set('');
    this.diversityData = [];
    this.availableClusters = [];
    this.selectedClusters = [];
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

    console.log('Starting cluster diversity analysis with parameters:', params);

    // Call the cluster diversity API
    this.apiService.clusterDiversity(params).pipe(
      tap(response => {
        if (response.status === 'ok' && response.result) {
          this.processedFileName.set(response.result);
          console.log('Cluster diversity analysis complete:', response.result);
        }
      }),
      switchMap(response => {
        // Automatically load results after successful analysis
        if (response.status === 'ok' && response.result) {
          return this.apiService.downloadFile(response.result).pipe(
            tap(blob => this.parseFileBlob(blob, response.result))
          );
        }
        return of(null);
      }),
      catchError(error => {
        const errorMsg = error.error?.detail || 'Cluster diversity analysis failed';
        console.error('Cluster diversity error:', errorMsg);
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
      this.parseResultFile(text, filename);
    };
    reader.readAsText(blob);
  }

  parseResultFile(content: string, filename: string): void {
    const isCsv = filename.endsWith('.csv');
    this.diversityData = [];
    
    if (isCsv) {
      // Parse CSV content
      const lines = content.split('\n').filter(line => line.trim());
      if (lines.length > 1) {
        const headers = lines[0].split(',');
        
        for (let i = 1; i < lines.length; i++) {
          const values = lines[i].split(',');
          if (values.length >= 7) {
            const clusterData = {
              cluster: parseInt(values[0]) || 0,
              totalSequences: parseInt(values[1]) || 0,
              totalReads: parseInt(values[2]) || 0,
              totalRPU: parseFloat(values[3]) || 0,
              averageLED: parseFloat(values[4]) || 0,
              SID: parseFloat(values[5]) || 0,
              seedID: values[6] || '',
              seedSequence: values[7] || ''
            };
            this.diversityData.push(clusterData);
          }
        }
      }
    }

    // Update available clusters for k-mer plot selection
    this.availableClusters = this.diversityData.map(d => d.cluster).sort((a, b) => a - b);
    
    // Emit results to parent component for display in right section
    this.resultsReady.emit({
      data: this.diversityData,
      inputFile: this.savedFileName
    });
    
    console.log('Parsed diversity data:', this.diversityData.length, 'clusters');
  }

  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) {
      console.warn('No file available for download. Please run cluster diversity analysis first.');
      return;
    }

    console.log('Downloading file:', filename);

    this.apiService.downloadFile(filename).subscribe({
      next: (blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        window.URL.revokeObjectURL(url);
        console.log('Download started');
      },
      error: (error) => {
        const errorMsg = error.error?.detail || 'Download failed';
        console.error('Download error:', errorMsg);
      }
    });
  }

  clusterMetaPlot(): void {
    if (this.diversityData.length === 0) {
      alert('No data available for plotting. Please run the diversity analysis first.');
      return;
    }

    const params = {
      xaxis: this.clusterMetaXAxis,
      yaxes: [this.clusterMetaYAxis1, this.clusterMetaYAxis2, this.clusterMetaYAxis3],
      plot_title: this.clusterMetaPlotTitle,
      line_colours: [this.clusterMetaColor1, this.clusterMetaColor2, this.clusterMetaColor3]
    };

    this.showMetaplotModal.emit({ data: this.diversityData, params });
  }

  kmerPlot(): void {
    if (this.diversityData.length === 0) {
      alert('No data available for plotting. Please run the diversity analysis first.');
      return;
    }

    if (this.selectedClusters.length === 0) {
      alert('Please select at least one cluster for the k-mer plot.');
      return;
    }

    // Prepare k-mer analysis request
    const kmerRequest: KmerAnalysisRequest = {
      input_path: this.savedFileName,
      k_size: parseInt(this.kmerSize),
      selected_clusters: this.selectedClusters,
      method: this.plotType.toLowerCase() as 'pca' | 'umap'
    };

    console.log('Starting k-mer analysis with parameters:', kmerRequest);

    // Prepare parameters for the plot modal
    const params = {
      kmerSize: parseInt(this.kmerSize),
      clustersToPlot: this.selectedClusters,
      plotType: this.plotType,
      xaxis: this.kmerPlotXAxis,
      yaxis: this.kmerPlotYAxis,
      legend_title: this.kmerPlotLegendTitle,
      plot_title: this.kmerPlotTitle,
      colour_palette_disc: this.kmerPlotColorPalette,
      kmerRequest: kmerRequest
    };

    this.showKmerPlotModal.emit({ data: this.diversityData, params });
  }

  // Helper methods for improved cluster selection UI
  getAvailableClustersForSelection(): number[] {
    return this.availableClusters.filter(cluster => !this.selectedClusters.includes(cluster));
  }

  addSelectedCluster(cluster: number): void {
    if (!this.selectedClusters.includes(cluster)) {
      this.selectedClusters.push(cluster);
      this.selectedClusters.sort((a, b) => a - b); // Keep sorted
    }
  }

  removeSelectedCluster(cluster: number): void {
    this.selectedClusters = this.selectedClusters.filter(c => c !== cluster);
  }

  // Helper method to track cluster selection
  onClusterSelectionChange(clusters: number[]): void {
    this.selectedClusters = clusters;
  }
}
