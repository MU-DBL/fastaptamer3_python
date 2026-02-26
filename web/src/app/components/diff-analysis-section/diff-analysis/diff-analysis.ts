import { Component, inject, signal, ChangeDetectorRef, NgZone } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { ApiService } from '../../../shared/api.service';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { Table, TableConfig } from '../../common/table/table';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';

interface UploadingFile {
  name: string;
  progress: number;
  isComplete: boolean;
}

@Component({
  selector: 'app-diff-analysis',
  imports: [
    CommonModule,
    FormsModule,
    Table,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './diff-analysis.html',
  styleUrl: './diff-analysis.scss',
  standalone: true
})
export class DiffAnalysis {
  // Services
  private apiService = inject(ApiService);
  private plotModalService = inject(PlotModalService);
  private cdr = inject(ChangeDetectorRef);
  private ngZone = inject(NgZone);

  // Signals
  isProcessing = signal(false);
  processedFileName = signal('');

  // File handling - Condition 1
  uploadedFilesCond1: string[] = [];
  uploadingFilesCond1: UploadingFile[] = [];
  isUploadingCond1: boolean = false;

  // File handling - Condition 2
  uploadedFilesCond2: string[] = [];
  uploadingFilesCond2: UploadingFile[] = [];
  isUploadingCond2: boolean = false;

  // Differential analysis parameters
  pCutoff: number = 0.1;
  downloadFormat: string = 'csv';

  // Data State
  diffAnalysisData: any[] = [];

  // Table configuration
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

  // edgeR plot customization
  adjustEdgeRPlot: string = 'no';
  edgeRPlotXAxis: string = 'logCPM';
  edgeRPlotYAxis: string = 'logFC';
  edgeRPlotTitle: string = 'edgeR results';
  edgeRPlotSigColor: string = '#FF0000';
  edgeRPlotInsigColor: string = '#000000';

  // Watch for changes in adjustment toggle
  onAdjustEdgeRPlotChange(): void {
    if (this.adjustEdgeRPlot === 'no') {
      this.edgeRPlotXAxis = 'logCPM';
      this.edgeRPlotYAxis = 'logFC';
      this.edgeRPlotTitle = 'edgeR results';
      this.edgeRPlotSigColor = '#FF0000';
      this.edgeRPlotInsigColor = '#000000';
    }
  }

  // Helper method for slider label formatting
  formatPValue(value: number): string {
    return value.toFixed(2);
  }

  // ========================================================================
  // FILE UPLOAD HANDLERS - CONDITION 1
  // ========================================================================
  
  onCondition1FilesSelected(event: any): void {
    const files: FileList = event.target.files;
    if (!files || files.length === 0) return;

    // Clear previous uploads
    this.uploadedFilesCond1 = [];
    this.uploadingFilesCond1 = [];

    // Initialize upload tracking
    this.uploadingFilesCond1 = Array.from(files).map(file => ({
      name: file.name,
      progress: 0,
      isComplete: false
    }));

    this.isUploadingCond1 = true;
    console.log('Condition 1 files selected:', files.length);
    
    // Upload files sequentially
    this.uploadFilesSequentially(Array.from(files), 0, 'cond1');
  }

  // ========================================================================
  // FILE UPLOAD HANDLERS - CONDITION 2
  // ========================================================================
  
  onCondition2FilesSelected(event: any): void {
    const files: FileList = event.target.files;
    if (!files || files.length === 0) return;

    // Clear previous uploads
    this.uploadedFilesCond2 = [];
    this.uploadingFilesCond2 = [];

    // Initialize upload tracking
    this.uploadingFilesCond2 = Array.from(files).map(file => ({
      name: file.name,
      progress: 0,
      isComplete: false
    }));

    this.isUploadingCond2 = true;
    console.log('Condition 2 files selected:', files.length);
    
    // Upload files sequentially
    this.uploadFilesSequentially(Array.from(files), 0, 'cond2');
  }

  // ========================================================================
  // SEQUENTIAL FILE UPLOAD
  // ========================================================================
  
  uploadFilesSequentially(files: File[], index: number, condition: 'cond1' | 'cond2'): void {
    if (index >= files.length) {
      // All files uploaded
      if (condition === 'cond1') {
        this.isUploadingCond1 = false;
      } else {
        this.isUploadingCond2 = false;
      }
      console.log(`All ${condition} files uploaded successfully`);
      return;
    }

    const file = files[index];
    const uploadingFiles = condition === 'cond1' ? this.uploadingFilesCond1 : this.uploadingFilesCond2;
    const uploadingFile = uploadingFiles[index];

    // Start progress simulation
    const progressInterval = setInterval(() => {
      if (uploadingFile.progress < 90) {
        this.ngZone.run(() => {
          uploadingFile.progress += 10;
          this.cdr.markForCheck();
        });
      }
    }, 100);

    // Upload to backend
    this.apiService.uploadFile(file).subscribe({
      next: (response) => {
        this.ngZone.run(() => {
          clearInterval(progressInterval);
          uploadingFile.progress = 100;
          uploadingFile.isComplete = true;
          
          // Add to uploaded files list
          if (condition === 'cond1') {
            this.uploadedFilesCond1.push(response.saved_filename);
          } else {
            this.uploadedFilesCond2.push(response.saved_filename);
          }
          
          console.log(`${condition} file ${index + 1} uploaded:`, response.saved_filename);
          
          // Trigger change detection
          this.cdr.markForCheck();
          
          // Continue with next file
          setTimeout(() => {
            this.uploadFilesSequentially(files, index + 1, condition);
          }, 200);
        });
      },
      error: (error) => {
        clearInterval(progressInterval);
        uploadingFile.progress = 0;
        uploadingFile.isComplete = false;
        console.error(`Failed to upload ${condition} file ${index + 1}:`, error);
        alert(`Failed to upload ${file.name}`);
        
        // Continue with next file despite error
        this.uploadFilesSequentially(files, index + 1, condition);
      }
    });
  }

  // ========================================================================
  // DIFFERENTIAL ANALYSIS
  // ========================================================================
  
  onStart(): void {
    if (this.uploadedFilesCond1.length < 2) {
      alert('Please upload at least 2 files for Condition 1.');
      return;
    }

    if (this.uploadedFilesCond2.length < 2) {
      alert('Please upload at least 2 files for Condition 2.');
      return;
    }

    this.isProcessing.set(true);
    this.diffAnalysisData = [];

    const params = {
      cond1_paths: this.uploadedFilesCond1,
      cond2_paths: this.uploadedFilesCond2,
      p_cutoff: this.pCutoff,
      output_format: this.downloadFormat
    };

    console.log('Starting differential analysis with parameters:', params);

    this.apiService.differentialAnalysis(params).pipe(
      tap(response => {
        console.log('Differential analysis response:', response);
        this.processedFileName.set(response.result);
      }),
      switchMap(response => {
        if (response.result) {
          return this.apiService.downloadFile(response.result).pipe(
            tap(blob => this.parseFileBlob(blob, response.result))
          );
        }
        throw new Error('No result file returned from differential analysis');
      }),
      catchError(error => {
        console.error('Differential analysis error:', error);
        const errorMsg = error.error?.detail || error.message || 'Unknown error';
        alert(`Differential analysis failed: ${errorMsg}`);
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
    
    this.diffAnalysisData = [];
    
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
          
          // Convert numeric fields
          if (['logFC', 'logCPM', 'PValue'].includes(header)) {
            row[header] = value ? parseFloat(value) : 0;
          } else {
            row[header] = value || '';
          }
        });
        
        this.diffAnalysisData.push(row);
      }
    }
    
    console.log('Parsed differential analysis data:', this.diffAnalysisData.length, 'rows');
    this.cdr.detectChanges();
  }

  onDownload(): void {
    if (!this.processedFileName()) {
      alert('No differential analysis data available to download.');
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
  // EDGER PLOT
  // ========================================================================
  
  async edgeRPlot(): Promise<void> {
    if (this.diffAnalysisData.length === 0) {
      alert('No differential analysis data available for plotting. Please run the analysis first.');
      return;
    }

    // Extract data for plotting
    const logCPM = this.diffAnalysisData.map(row => row['logCPM']);
    const logFC = this.diffAnalysisData.map(row => row['logFC']);
    const pClass = this.diffAnalysisData.map(row => row['PClass']);
    const sequences = this.diffAnalysisData.map(row => row['Sequence'] || '');

    // Separate significant and insignificant points
    const sigIndices = pClass.map((p, i) => p === 'Sig.' ? i : -1).filter(i => i >= 0);
    const insigIndices = pClass.map((p, i) => p === 'Not Sig.' ? i : -1).filter(i => i >= 0);

    const traces = [];

    // Insignificant points
    if (insigIndices.length > 0) {
      traces.push({
        x: insigIndices.map(i => logCPM[i]),
        y: insigIndices.map(i => logFC[i]),
        mode: 'markers',
        type: 'scatter',
        marker: {
          color: this.edgeRPlotInsigColor,
          size: 6,
          opacity: 0.5
        },
        text: insigIndices.map(i => sequences[i]),
        hovertemplate: '%{text}<br>logCPM: %{x:.3f}<br>logFC: %{y:.3f}<extra></extra>',
        name: 'Not Significant'
      });
    }

    // Significant points
    if (sigIndices.length > 0) {
      traces.push({
        x: sigIndices.map(i => logCPM[i]),
        y: sigIndices.map(i => logFC[i]),
        mode: 'markers',
        type: 'scatter',
        marker: {
          color: this.edgeRPlotSigColor,
          size: 6,
          opacity: 0.7
        },
        text: sigIndices.map(i => sequences[i]),
        hovertemplate: '%{text}<br>logCPM: %{x:.3f}<br>logFC: %{y:.3f}<extra></extra>',
        name: 'Significant'
      });
    }

    const layout = {
      title: {
        text: this.edgeRPlotTitle,
        font: { size: 18, family: 'Arial, sans-serif', weight: 'bold' }
      },
      xaxis: {
        title: {
          text: this.edgeRPlotXAxis,
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
          text: this.edgeRPlotYAxis,
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
      hovermode: 'closest',
      showlegend: true
    };

    this.plotModalService.openPlot({
      data: traces,
      layout: layout,
      config: { responsive: true }
    });
  }
}
