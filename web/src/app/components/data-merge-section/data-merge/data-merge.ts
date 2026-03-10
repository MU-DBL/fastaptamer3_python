import { Component, inject, signal, ChangeDetectorRef, OnDestroy } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { ApiService } from '../../../shared/api.service';
import { FileService } from '../../../shared/file-service';
import { PlotModalService } from '../../../shared/plot-modal.service';
import { Table, TableConfig } from '../../common/table/table';
import { CdkDragDrop, moveItemInArray, DragDropModule } from '@angular/cdk/drag-drop';
import { switchMap, tap, catchError, finalize } from 'rxjs/operators';
import { of } from 'rxjs';

interface FileSelection {
  fileName: string;
  order: number | null;
  include: boolean;
}

interface UploadingFile {
  name: string;
  progress: number;
  isComplete: boolean;
  savedFileName?: string;
}

@Component({
  selector: 'app-data-merge',
  imports: [
    CommonModule,
    FormsModule,
    Table,
    DragDropModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './data-merge.html',
  styleUrl: './data-merge.scss',
  standalone: true
})
export class DataMerge implements OnDestroy {
  // Services
  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private plotModalService = inject(PlotModalService);
  private cdr = inject(ChangeDetectorRef);

  // File management
  uploadedFiles: string[] = [];
  fileSelections: FileSelection[] = [];
  uploadingFiles: UploadingFile[] = [];
  isUploading: boolean = false;

  // Merge parameters
  mergeType: 'Union' | 'Intersection' | 'Left' = 'Union';

  // Processing state
  isProcessing = signal(false);
  isGeneratingPlot = signal(false);
  mergedFileName = signal('');
  mergedData: any[] = [];

  // Table configuration
  tableConfig: TableConfig = {
    columns: [],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  // Persistence plot customization
  adjustPersistencePlot: string = 'no';
  persistencePlotXAxis: string = 'Population count';
  persistencePlotYAxis: string = 'Unique read count';
  persistencePlotTitle: string = 'Sequence persistence analysis';
  persistenceBarOutline: string = '#000000';
  persistenceBarFill: string = '#87CEEB';

  // UpSet plot customization
  adjustUpSetPlot: string = 'no';
  upsetPlotXAxis: string = 'Sequences per set';
  upsetPlotYAxis: string = 'Sequence intersections';
  upsetBarFill: string = '#87CEEB';

  // ========================================================================
  // FILE UPLOAD AND MANAGEMENT
  // ========================================================================

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
  }

  ngOnDestroy(): void {
    this.uploadedFiles.forEach(f => this.apiService.deleteFile(f).subscribe());
    if (this.mergedFileName()) {
      this.apiService.deleteFile(this.mergedFileName()).subscribe();
    }
  }

  onMultipleFilesSelected(event: any): void {
    const files: FileList = event.target.files;
    if (!files || files.length === 0) return;

    // Clear previous upload progress
    this.uploadingFiles = [];
    this.mergedFileName.set('');
    this.mergedData = [];

    // Initialize upload tracking
    this.uploadingFiles = Array.from(files).map(file => ({
      name: file.name,
      progress: 0,
      isComplete: false
    }));

    this.isUploading = true;
    console.log('Multiple files selected:', files.length);

    // Upload files sequentially
    this.uploadFilesSequentially(Array.from(files), 0);
  }

  uploadFilesSequentially(files: File[], index: number): void {
    if (index >= files.length) {
      this.isUploading = false;
      console.log('All files uploaded successfully');
      
      setTimeout(() => {
        this.autoPopulateFileSelections();
        this.cdr.markForCheck();
      }, 0);
      return;
    }

    const file = files[index];
    const uploadingFile = this.uploadingFiles[index];

    // Simulate progress
    const progressInterval = setInterval(() => {
      if (uploadingFile.progress < 90) {
        uploadingFile.progress += 10;
      }
    }, 50);

    this.apiService.uploadFile(file).subscribe({
      next: (response) => {
        clearInterval(progressInterval);
        uploadingFile.progress = 100;
        uploadingFile.isComplete = true;
        uploadingFile.savedFileName = response.saved_filename;
        
        this.uploadedFiles.push(response.saved_filename);
        console.log(`File ${index + 1} uploaded:`, response.saved_filename);
        
        setTimeout(() => {
          this.uploadFilesSequentially(files, index + 1);
        }, 200);
      },
      error: (error) => {
        clearInterval(progressInterval);
        uploadingFile.progress = 0;
        uploadingFile.isComplete = false;
        
        console.error(`Failed to upload file ${index + 1}:`, error);
        alert(`Failed to upload ${file.name}`);
        
        this.uploadFilesSequentially(files, index + 1);
      }
    });
  }

  autoPopulateFileSelections(): void {
    this.fileSelections = this.uploadedFiles.map((fileName, index) => ({
      fileName: fileName,
      order: index + 1,
      include: true
    }));
    console.log('Auto-populated file selections:', this.fileSelections.length);
  }

  // Drag and drop handler
  onFileDrop(event: CdkDragDrop<FileSelection[]>): void {
    moveItemInArray(this.fileSelections, event.previousIndex, event.currentIndex);
    this.updateOrderNumbers();
  }

  moveFileUp(index: number): void {
    if (index > 0) {
      const temp = this.fileSelections[index];
      this.fileSelections[index] = this.fileSelections[index - 1];
      this.fileSelections[index - 1] = temp;
      this.updateOrderNumbers();
    }
  }

  moveFileDown(index: number): void {
    if (index < this.fileSelections.length - 1) {
      const temp = this.fileSelections[index];
      this.fileSelections[index] = this.fileSelections[index + 1];
      this.fileSelections[index + 1] = temp;
      this.updateOrderNumbers();
    }
  }

  updateOrderNumbers(): void {
    this.fileSelections.forEach((selection, index) => {
      selection.order = index + 1;
    });
  }

  getOrderedFiles(): string[] {
    return this.fileSelections
      .filter(fs => fs.fileName && fs.order !== null && fs.include)
      .sort((a, b) => (a.order || 0) - (b.order || 0))
      .map(fs => fs.fileName);
  }

  // ========================================================================
  // MERGE PROCESSING
  // ========================================================================

  validateInputs(): boolean {
    const orderedFiles = this.getOrderedFiles();
    
    if (orderedFiles.length < 2) {
      alert('Please include at least 2 files for merging!');
      return false;
    }

    return true;
  }

  onStart(): void {
    if (!this.validateInputs()) {
      return;
    }

    this.isProcessing.set(true);
    this.mergedFileName.set('');
    this.mergedData = [];

    const orderedFiles = this.getOrderedFiles();
    
    // Map UI merge type to backend merge type
    const mergeTypeMap: { [key: string]: string } = {
      'Union': 'outer',
      'Intersection': 'inner',
      'Left': 'left'
    };

    const params = {
      input_paths: orderedFiles,
      merge_type: mergeTypeMap[this.mergeType],
      output_format: 'csv'
    };

    this.apiService.dataMerge(params).pipe(
      tap(response => {
        if (response.status === 'ok' && response.result) {
          this.mergedFileName.set(response.result);
          console.log('Merge complete:', response.result);
        }
      }),
      switchMap(response => {
        if (response.status === 'ok' && response.result) {
          return this.apiService.downloadFile(response.result).pipe(
            switchMap(blob => this.fileService.parseClusterFile(blob, response.result)),
            tap(parsedData => {
              this.mergedData = parsedData;
              this.updateTableConfig(parsedData);
            })
          );
        }
        return of(null);
      }),
      catchError(error => {
        const errorMsg = error.error?.detail || 'Merge failed';
        console.error('Merge error:', errorMsg);
        alert(`Error: ${errorMsg}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  updateTableConfig(data: any[]): void {
    if (data.length === 0) return;

    const keys = Object.keys(data[0]);
    this.tableConfig = {
      ...this.tableConfig,
      columns: keys.map(key => ({
        key: key,
        label: key
      }))
    };
  }

  // ========================================================================
  // DOWNLOAD
  // ========================================================================

  onDownload(): void {
    const filename = this.mergedFileName();
    if (!filename) {
      console.warn('No file available for download.');
      return;
    }
    this.fileService.downloadFile(filename);
  }

  // ========================================================================
  // PLOT CUSTOMIZATIONS
  // ========================================================================

  onAdjustPersistencePlotChange(): void {
    if (this.adjustPersistencePlot === 'no') {
      this.resetPersistencePlotDefaults();
    }
  }

  resetPersistencePlotDefaults(): void {
    this.persistencePlotXAxis = 'Population count';
    this.persistencePlotYAxis = 'Unique read count';
    this.persistencePlotTitle = 'Sequence persistence analysis';
    this.persistenceBarOutline = '#000000';
    this.persistenceBarFill = '#87CEEB';
  }

  onAdjustUpSetPlotChange(): void {
    if (this.adjustUpSetPlot === 'no') {
      this.resetUpSetPlotDefaults();
    }
  }

  resetUpSetPlotDefaults(): void {
    this.upsetPlotXAxis = 'Sequences per set';
    this.upsetPlotYAxis = 'Sequence intersections';
    this.upsetBarFill = '#87CEEB';
  }

  // ========================================================================
  // PERSISTENCE PLOT
  // ========================================================================

  async persistencePlot(): Promise<void> {
    if (!this.mergedFileName()) {
      alert('Please merge files first before generating plots.');
      return;
    }

    this.isGeneratingPlot.set(true);

    const params = {
      merged_file_path: this.mergedFileName()
    };

    this.apiService.sequencePersistence(params).subscribe({
      next: (response) => {
        if (response.status === 'ok' && response.data) {
          this.createPersistencePlot(response.data);
        }
      },
      error: (error) => {
        const errorMsg = error.error?.detail || 'Failed to generate persistence plot';
        console.error('Persistence plot error:', errorMsg);
        alert(`Error: ${errorMsg}`);
      },
      complete: () => {
        this.isGeneratingPlot.set(false);
      }
    });
  }

  createPersistencePlot(data: Array<{ freq: number; seqCount: number }>): void {
    const xData = data.map(d => d.freq);
    const yData = data.map(d => d.seqCount);

    const trace = {
      x: xData,
      y: yData,
      type: 'bar',
      marker: {
        color: this.persistenceBarFill,
        line: {
          color: this.persistenceBarOutline,
          width: 1
        }
      },
      text: yData.map(String),
      textposition: 'auto',
      hovertemplate: '<b>Population Count:</b> %{x}<br><b>Unique Sequences:</b> %{y}<extra></extra>'
    };

    const layout = {
      title: {
        text: this.persistencePlotTitle,
        font: { size: 18 }
      },
      xaxis: {
        title: {
          text: this.persistencePlotXAxis,
          font: { size: 14 }
        },
        showgrid: true,
        gridcolor: '#e0e0e0',
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        mirror: true,
        tickmode: 'linear',
        tick0: 0,
        dtick: 1
      },
      yaxis: {
        title: {
          text: this.persistencePlotYAxis,
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
      margin: { t: 80, b: 90, l: 100, r: 60 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white',
      bargap: 0.2
    };

    this.plotModalService.openPlot({
      data: [trace],
      layout: layout,
      config: { responsive: true }
    });
  }

  // ========================================================================
  // UPSET PLOT
  // ========================================================================

  async upsetPlot(): Promise<void> {
    if (!this.mergedFileName()) {
      alert('Please merge files first before generating plots.');
      return;
    }

    this.isGeneratingPlot.set(true);

    // Get original file names for set labels
    const orderedFiles = this.getOrderedFiles();
    const setNames = orderedFiles.map(filename => {
      // Remove file extension and path for cleaner labels
      const name = filename.split('/').pop() || filename;
      return name.replace(/\.(fasta|fa|csv|tsv)$/i, '');
    });

    const params = {
      merged_file_path: this.mergedFileName(),
      fasta_names: setNames
    };

    this.apiService.upSetData(params).subscribe({
      next: (response) => {
        if (response.status === 'ok') {
          this.createUpSetPlot(response);
        }
      },
      error: (error) => {
        const errorMsg = error.error?.detail || 'Failed to generate UpSet plot';
        console.error('UpSet plot error:', errorMsg);
        alert(`Error: ${errorMsg}`);
      },
      complete: () => {
        this.isGeneratingPlot.set(false);
      }
    });
  }

  createUpSetPlot(data: {
    sets: string[];
    set_sizes: { [key: string]: number };
    intersections: Array<{ sets: string[]; size: number; sequences: string[] }>;
    total_unique_sequences: number;
  }): void {
    // Create enhanced bar chart showing intersection sizes with color coding
    this.createUpSetBarChart(data);
  }

  private createUpSetBarChart(data: {
    sets: string[];
    set_sizes: { [key: string]: number };
    intersections: Array<{ sets: string[]; size: number; sequences: string[] }>;
    total_unique_sequences: number;
  }): void {
    // Enhanced bar chart showing intersection sizes
    const sortedIntersections = data.intersections
      .sort((a, b) => b.size - a.size)
      .slice(0, 30); // Show top 30 intersections

    // Create labels showing set combinations with better formatting
    const xLabels = sortedIntersections.map(i => {
      if (i.sets.length === 1) {
        return `[${i.sets[0]}]`;
      }
      return i.sets.join(' ∩ ');
    });
    const yData = sortedIntersections.map(i => i.size);

    // Color bars by number of sets in intersection
    const colors = sortedIntersections.map(i => {
      const numSets = i.sets.length;
      if (numSets === 1) return '#2196f3'; // Single set - blue
      if (numSets === 2) return '#4caf50'; // Two sets - green  
      if (numSets === 3) return '#ff9800'; // Three sets - orange
      return '#f44336'; // More sets - red
    });

    const trace = {
      x: xLabels,
      y: yData,
      type: 'bar',
      marker: {
        color: colors,
        line: {
          color: '#000',
          width: 1
        }
      },
      text: yData.map(String),
      textposition: 'outside',
      hovertemplate: '<b>Sets:</b> %{x}<br><b>Sequences:</b> %{y}<extra></extra>'
    };

    const layout = {
      title: {
        text: 'UpSet Plot - Top 30 Intersections',
        font: { size: 18, weight: 600 }
      },
      xaxis: {
        title: {
          text: this.upsetPlotYAxis,
          font: { size: 14, weight: 600 }
        },
        showgrid: true,
        gridcolor: '#e0e0e0',
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        mirror: true,
        tickangle: -45,
        automargin: true
      },
      yaxis: {
        title: {
          text: this.upsetPlotXAxis,
          font: { size: 14, weight: 600 }
        },
        showgrid: true,
        gridcolor: '#e0e0e0',
        showline: true,
        linewidth: 2,
        linecolor: 'black',
        mirror: true,
        automargin: true
      },
      autosize: true,
      height: 650,
      margin: { t: 80, b: 180, l: 100, r: 60 },
      plot_bgcolor: 'white',
      paper_bgcolor: 'white',
      showlegend: false
    };

    this.plotModalService.openPlot({
      data: [trace],
      layout: layout,
      config: { responsive: true, displayModeBar: true }
    });
  }
}
