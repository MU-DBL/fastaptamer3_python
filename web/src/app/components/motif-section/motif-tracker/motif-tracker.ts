import { Component, inject, signal, OnDestroy } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { switchMap, tap, catchError, finalize, map } from 'rxjs/operators';
import { of } from 'rxjs';
import { Table, TableConfig } from '../../common/table/table';
import { CdkDragDrop, moveItemInArray, DragDropModule } from '@angular/cdk/drag-drop';
import { TrackerPlot } from './charts/tracker-plot';
import { SplitPanel } from '../../common/split-panel/split-panel';

interface FileSelection {
  fileName: string;
  order: number | null;
  include: boolean; // Whether to include this file in tracking
}

@Component({
  selector: 'app-motif-tracker',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    Table,
    TrackerPlot,
    SplitPanel,
    DragDropModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './motif-tracker.html',
  styleUrl: './motif-tracker.scss'
})
export class MotifTracker implements OnDestroy {
  private apiService = inject(ApiService);

  // Table configurations
  motifTrackerTableConfig: TableConfig = {
    columns: [
      { key: 'PopulationNumber', label: 'Pop #' },
      { key: 'PopulationName', label: 'Population' },
      { key: 'Motif', label: 'Motif' },
      { key: 'Alias', label: 'Alias' },
      { key: 'TotalReads', label: 'Total Reads' },
      { key: 'Percentage', label: 'Percentage' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  sequenceTrackerTableConfig: TableConfig = {
    columns: [
      { key: 'PopulationNumber', label: 'Pop #' },
      { key: 'PopulationName', label: 'Population' },
      { key: 'sequences', label: 'Sequence' },
      { key: 'Alias', label: 'Alias' },
      { key: 'Rank', label: 'Rank' },
      { key: 'Reads', label: 'Reads' },
      { key: 'RPU', label: 'RPU' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  enrichmentTableConfig: TableConfig = {
    columns: [
      { key: 'Comparison', label: 'Comparison' },
      { key: 'Query', label: 'Query' },
      { key: 'Alias', label: 'Alias' },
      { key: 'Enrichment', label: 'Enrichment' }
    ],
    initialPageSize: 10,
    pageSizeOptions: [10, 25, 50, 100]
  };

  // Get dynamic table config based on query type
  get trackerTableConfig(): TableConfig {
    return this.queryType === 'Motif' ? this.motifTrackerTableConfig : this.sequenceTrackerTableConfig;
  }

  // Plot state and parameters
  showTrackerPlot = signal(false);
  
  // Plot customization values
  plotXAxis: string = 'Population';
  plotYAxis: string = 'Percentage';
  plotTitle: string = 'Motif tracker';
  plotColorPalette: string = 'Dark2';
  adjustPlotSettings: string = 'no';

  // Get default values based on query type
  get defaultPlotYAxis(): string {
    return this.queryType === 'Motif' ? 'Percentage' : 'RPU';
  }

  get defaultPlotTitle(): string {
    return this.queryType === 'Motif' ? 'Motif tracker' : 'Sequence tracker';
  }

  // File management
  uploadedFiles: string[] = [];
  fileSelections: FileSelection[] = [];
  
  // Form values
  queryList: string = '';
  queryAliases: string = '';
  queryType: 'Motif' | 'Sequence' = 'Motif';
  motifType: string = 'Nucleotide'; // Only for Motif tracking

  // Processing states
  isProcessing = signal(false);
  processedFileName = signal('');
  enrichmentFileName = signal('');

  // Table data
  trackerData: any[] = [];
  enrichmentData: any[] = [];

  cancelProcessing(): void {
    this.apiService.cancelProcesses().subscribe();
    this.isProcessing.set(false);
  }

  ngOnDestroy(): void {
    this.uploadedFiles.forEach(f => this.apiService.deleteFile(f).subscribe());
    if (this.processedFileName()) {
      this.apiService.deleteFile(this.processedFileName()).subscribe();
    }
    if (this.enrichmentFileName()) {
      this.apiService.deleteFile(this.enrichmentFileName()).subscribe();
    }
  }

  // ========================================================================
  // FILE MANAGEMENT
  // ========================================================================

  onUploadComplete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadedFiles.push(result.savedFileName);
      this.fileSelections.push({
        fileName: result.savedFileName,
        order: this.fileSelections.length + 1,
        include: true,
      });
    }
  }

  removeFile(selection: FileSelection): void {
    this.apiService.deleteFile(selection.fileName).subscribe();
    this.uploadedFiles = this.uploadedFiles.filter(f => f !== selection.fileName);
    this.fileSelections = this.fileSelections.filter(s => s !== selection);
    this.updateOrderNumbers();
  }

  // Drag and drop handler
  onFileDrop(event: CdkDragDrop<FileSelection[]>): void {
    moveItemInArray(this.fileSelections, event.previousIndex, event.currentIndex);
    this.updateOrderNumbers();
  }

  moveFileUp(index: number): void {
    if (index > 0) {
      // Swap with previous item
      const temp = this.fileSelections[index];
      this.fileSelections[index] = this.fileSelections[index - 1];
      this.fileSelections[index - 1] = temp;
      
      // Update order numbers
      this.updateOrderNumbers();
    }
  }

  moveFileDown(index: number): void {
    if (index < this.fileSelections.length - 1) {
      // Swap with next item
      const temp = this.fileSelections[index];
      this.fileSelections[index] = this.fileSelections[index + 1];
      this.fileSelections[index + 1] = temp;
      
      // Update order numbers
      this.updateOrderNumbers();
    }
  }

  removeFileSelection(index: number): void {
    this.fileSelections.splice(index, 1);
    this.updateOrderNumbers();
  }

  updateOrderNumbers(): void {
    // Reassign order numbers sequentially
    this.fileSelections.forEach((selection, index) => {
      selection.order = index + 1;
    });
  }

  getOrderedFiles(): string[] {
    return this.fileSelections
      .filter(fs => fs.fileName && fs.order !== null && fs.include) // Only include files marked as included
      .sort((a, b) => (a.order || 0) - (b.order || 0))
      .map(fs => fs.fileName);
  }

  // ========================================================================
  // FORM VALIDATION & PROCESSING
  // ========================================================================

  validateInputs(): boolean {
    const orderedFiles = this.getOrderedFiles();
    
    if (orderedFiles.length < 2) {
      alert('Please include at least 2 files for tracking!');
      return false;
    }

    // Check for duplicate orders among included files
    const includedSelections = this.fileSelections.filter(fs => fs.include && fs.order !== null);
    const orders = includedSelections.map(fs => fs.order);
    const uniqueOrders = new Set(orders);
    if (orders.length !== uniqueOrders.size) {
      alert('Please ensure each included file has a unique order number!');
      return false;
    }

    if (!this.queryList.trim()) {
      alert('Please enter at least one query motif or sequence!');
      return false;
    }

    // Validate query aliases if provided
    const queries = this.queryList.split('\n').filter(q => q.trim());
    const aliases = this.queryAliases.split('\n').filter(a => a.trim());
    
    if (aliases.length > 0 && aliases.length !== queries.length) {
      alert('When aliases are provided, there must be one for each query!');
      return false;
    }

    return true;
  }

  onStart(): void {
    if (!this.validateInputs()) {
      return;
    }

    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.enrichmentFileName.set('');
    this.trackerData = [];
    this.enrichmentData = [];

    const orderedFiles = this.getOrderedFiles();
    const queries = this.queryList.split('\n').filter(q => q.trim()).map(q => q.trim());
    const aliases = this.queryAliases.split('\n').filter(a => a.trim()).map(a => a.trim());

    // Use filenames as population names
    const populationNames = orderedFiles.map(f => f);

    const baseParams = {
      input_paths: orderedFiles,
      population_names: populationNames,
      query_list: queries,
      query_aliases: aliases.length > 0 ? aliases : undefined
    };

    const apiCall = this.queryType === 'Motif'
      ? this.apiService.motifTracker({ ...baseParams, motif_type: this.motifType })
      : this.apiService.sequenceTracker(baseParams);

    apiCall.pipe(
      tap(response => {
        if (response.status === 'ok') {
          this.processedFileName.set(response.result);
          this.enrichmentFileName.set(response.enrichment);
          console.log('Tracking complete:', response);
        }
      }),
      switchMap(response => {
        // Load tracker results
        if (response.status === 'ok' && response.result) {
          return this.apiService.fetchFileText(response.result).pipe(
            tap(text => this.parseTrackerContent(text, response.result))
          );
        }
        return of(null);
      }),
      switchMap(() => {
        // Load enrichment results
        if (this.enrichmentFileName()) {
          return this.apiService.fetchFileText(this.enrichmentFileName()).pipe(
            tap(text => this.parseEnrichmentContent(text))
          );
        }
        return of(null);
      }),
      catchError(error => {
        const errorMsg = error.error?.detail || 'Tracking failed';
        console.error('Tracking error:', errorMsg);
        alert(`Error: ${errorMsg}`);
        return of(null);
      }),
      finalize(() => {
        this.isProcessing.set(false);
      })
    ).subscribe();
  }

  // ========================================================================
  // FILE PARSING
  // ========================================================================

  parseTrackerContent(content: string, filename: string): void {
    const isCsv = filename.endsWith('.csv');
    this.trackerData = [];

    if (isCsv) {
      const lines = content.split('\n').filter(line => line.trim());
      if (lines.length < 2) return;

      const headers = lines[0].split(',').map(h => h.trim());
      
      for (let i = 1; i < lines.length; i++) {
        const values = lines[i].split(',').map(v => v.trim());
        const row: any = {};
        
        headers.forEach((header, idx) => {
          row[header] = values[idx] || '';
        });
        
        this.trackerData.push(row);
      }
    } else {
      // Parse FASTA format if needed
      console.log('FASTA format not implemented for tracker results');
    }

    console.log('Parsed tracker data:', this.trackerData.length, 'rows');
  }

  parseEnrichmentContent(content: string): void {
    this.enrichmentData = [];
    const lines = content.split('\n').filter(line => line.trim());
    if (lines.length < 2) return;

    const headers = lines[0].split(',').map(h => h.trim());
    
    for (let i = 1; i < lines.length; i++) {
      const values = lines[i].split(',').map(v => v.trim());
      const row: any = {};
      
      headers.forEach((header, idx) => {
        row[header] = values[idx] || '';
      });
      
      this.enrichmentData.push(row);
    }

    console.log('Parsed enrichment data:', this.enrichmentData.length, 'rows');
  }

  // ========================================================================
  // DOWNLOAD HANDLERS
  // ========================================================================

  onDownloadTracker(): void {
    if (!this.processedFileName()) {
      console.warn('No tracker file available for download');
      return;
    }

    const link = document.createElement('a');
    link.href = this.apiService.getDownloadUrl(this.processedFileName());
    link.download = this.processedFileName();
    link.click();
  }

  onDownloadEnrichment(): void {
    if (!this.enrichmentFileName()) {
      console.warn('No enrichment file available for download');
      return;
    }

    const link = document.createElement('a');
    link.href = this.apiService.getDownloadUrl(this.enrichmentFileName());
    link.download = this.enrichmentFileName();
    link.click();
  }

  // ========================================================================
  // PLOT METHODS
  // ========================================================================

  onAdjustPlotSettingsChange(): void {
    if (this.adjustPlotSettings === 'no') {
      // Reset to defaults when toggling back to "No"
      this.resetPlotSettings();
    }
  }

  resetPlotSettings(): void {
    this.plotXAxis = 'Population';
    this.plotYAxis = this.defaultPlotYAxis;
    this.plotTitle = this.defaultPlotTitle;
    this.plotColorPalette = 'Dark2';
  }

  // Update plot settings when query type changes
  onQueryTypeChange(): void {
    if (this.adjustPlotSettings === 'no') {
      // Auto-update defaults if user hasn't customized
      this.resetPlotSettings();
    }
  }

  openTrackerPlot(): void {
    if (this.trackerData.length === 0) {
      alert('No tracker results available for plotting. Please run the tracking process first.');
      return;
    }

    this.showTrackerPlot.set(true);
  }

  closeTrackerPlot(): void {
    this.showTrackerPlot.set(false);
  }
}

