import { Component, PLATFORM_ID, Inject, signal, WritableSignal } from '@angular/core';
import { isPlatformBrowser } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { ClusterDiversity } from '../cluster-diversity/cluster-diversity';
import { ClusterMsa } from '../cluster-msa/cluster-msa';
import { ClusterPhmm } from '../cluster-phmm/cluster-phmm';
import { Cluster } from '../cluster/cluster';
import { ApiService } from '../../../shared/api.service';
import { ClusterColorService } from '../../../shared/cluster-color.service';
import { 
  KmerAnalysisData, 
  KmerAnalysisRequest, 
  PlotlyTrace, 
  PlotlyLayout 
} from '../../../shared/kmer-analysis.types';
import { DiversityResultsEvent } from '../cluster-diversity/cluster-diversity';
import { ClusterResultsEvent } from '../cluster/cluster';
import { Recluster } from '../recluster/recluster';
import { PositionEnrichment } from '../position-enrichment/position-enrichment';

@Component({
  selector: 'app-clusterpage',
  imports: [
    CommonModule,
    FormsModule,
    Cluster,
    ClusterDiversity,
    ClusterMsa,
    ClusterPhmm,
    Recluster,
    PositionEnrichment,
    ...MATERIAL_IMPORTS
],
  templateUrl: './clusterpage.html',
  styleUrl: './clusterpage.scss',
})
export class Clusterpage {
  
  // Modal states
  showMetaplotModalFlag: boolean = false;
  showKmerPlotModalFlag: boolean = false;
  
  // Plot data and parameters
  metaplotData: any = null;
  metaplotParams: any = null;
  kmerPlotData: any = null;
  kmerPlotParams: any = null;
  
  // Active tab tracking (writable signal for two-way binding)
  selectedTabIndex: WritableSignal<number> = signal(0);
  
  // Tab-specific results storage
  tabResults: Map<number, {
    showResults: boolean;
    tableData: any[];
    filteredData: any[];
    paginatedData: any[];
    displayedColumns: string[];
    currentPage: number;
    pageSize: number;
    totalPages: number;
    searchText: string;
    columnFilters: any;
  }> = new Map();

  // Current tab results (computed from tabResults based on selectedTabIndex)
  get showResultsTable(): boolean {
    return this.tabResults.get(this.selectedTabIndex())?.showResults ?? false;
  }

  get tableData(): any[] {
    return this.tabResults.get(this.selectedTabIndex())?.tableData ?? [];
  }

  get filteredData(): any[] {
    return this.tabResults.get(this.selectedTabIndex())?.filteredData ?? [];
  }

  get paginatedData(): any[] {
    return this.tabResults.get(this.selectedTabIndex())?.paginatedData ?? [];
  }

  get displayedColumns(): string[] {
    return this.tabResults.get(this.selectedTabIndex())?.displayedColumns ?? 
           ['cluster', 'totalSequences', 'totalReads', 'totalRPU', 'averageLED', 'SID', 'seedID', 'seedSequence'];
  }

  get currentPage(): number {
    return this.tabResults.get(this.selectedTabIndex())?.currentPage ?? 1;
  }

  set currentPage(value: number) {
    const tabData = this.tabResults.get(this.selectedTabIndex());
    if (tabData) {
      tabData.currentPage = value;
    }
  }

  get pageSize(): number {
    return this.tabResults.get(this.selectedTabIndex())?.pageSize ?? 10;
  }

  set pageSize(value: number) {
    const tabData = this.tabResults.get(this.selectedTabIndex());
    if (tabData) {
      tabData.pageSize = value;
    }
  }

  get totalPages(): number {
    return this.tabResults.get(this.selectedTabIndex())?.totalPages ?? 0;
  }

  set totalPages(value: number) {
    const tabData = this.tabResults.get(this.selectedTabIndex());
    if (tabData) {
      tabData.totalPages = value;
    }
  }

  get searchText(): string {
    return this.tabResults.get(this.selectedTabIndex())?.searchText ?? '';
  }

  set searchText(value: string) {
    const tabData = this.tabResults.get(this.selectedTabIndex());
    if (tabData) {
      tabData.searchText = value;
    }
  }

  get columnFilters(): any {
    return this.tabResults.get(this.selectedTabIndex())?.columnFilters ?? {
      cluster: '',
      totalSequences: '',
      totalReads: '',
      totalRPU: '',
      averageLED: '',
      SID: '',
      seedID: '',
      seedSequence: '',
      sequenceID: '',  
      rank: '',           
      readCount: '',      
      rpu: '',            
      rankInCluster: '',  
      led: '',            
      sequence: ''
    };
  }

  set columnFilters(value: any) {
    const tabData = this.tabResults.get(this.selectedTabIndex());
    if (tabData) {
      tabData.columnFilters = value;
    }
  }
  
  // Plotly instance
  private plotly: any = null;
  
  // Input file tracking for k-mer analysis
  private lastProcessedFile: string = '';

  constructor(
    @Inject(PLATFORM_ID) private platformId: Object,
    private apiService: ApiService
  ) {
    // Load Plotly only on browser side
    if (isPlatformBrowser(this.platformId)) {
      this.loadPlotly();
    }
  }

  /**
   * Initialize or get tab-specific results
   */
  private ensureTabData(tabIndex: number): void {
    if (!this.tabResults.has(tabIndex)) {
      this.tabResults.set(tabIndex, {
        showResults: false,
        tableData: [],
        filteredData: [],
        paginatedData: [],
        displayedColumns: ['cluster', 'totalSequences', 'totalReads', 'totalRPU', 'averageLED', 'SID', 'seedID', 'seedSequence'],
        currentPage: 1,
        pageSize: 10,
        totalPages: 0,
        searchText: '',
        columnFilters: {
          cluster: '',
          totalSequences: '',
          totalReads: '',
          totalRPU: '',
          averageLED: '',
          SID: '',
          seedID: '',
          seedSequence: ''
        }
      });
    }
  }

  /**
   * Get current tab's data
   */
  private getCurrentTabData() {
    this.ensureTabData(this.selectedTabIndex());
    return this.tabResults.get(this.selectedTabIndex())!;
  }

  /**
   * Handle tab change event
   */
  onTabChange(index: number): void {
    this.selectedTabIndex.set(index);
  }

  /**
   * Handle diversity results from cluster-diversity component
   */
  onDiversityResultsReady(event: DiversityResultsEvent): void {
    const data = event.data;
    const inputFile = event.inputFile;
    
    const tabIndex = 1; // Diversity tab is index 1
    this.ensureTabData(tabIndex);
    const tabData = this.tabResults.get(tabIndex)!;
    
    tabData.tableData = data;
    tabData.filteredData = [...data];
    tabData.showResults = true;
    
    // Store the input filename for k-mer analysis
    if (inputFile) {
      this.lastProcessedFile = inputFile;
    }
    
    // Ensure we're on the correct tab to see results
    if (this.selectedTabIndex() !== tabIndex) {
      this.selectedTabIndex.set(tabIndex);
    }
    
    this.updatePagination();
    
    console.log('Diversity results ready:', data.length, 'clusters, showing on tab', tabIndex);
    console.log('Input file for k-mer analysis:', this.lastProcessedFile);
  }

  /**
   * Handle cluster results from cluster component
   */
  onClusterResultsReady(event: ClusterResultsEvent): void {
    const data = event.data;
    const inputFile = event.inputFile;
    
    const tabIndex = 0; // Cluster tab is index 0
    this.ensureTabData(tabIndex);
    const tabData = this.tabResults.get(tabIndex)!;
    
    // Update columns for cluster data structure
    tabData.displayedColumns = [
        'sequenceID', 
        'cluster', 
        'rankInCluster', 
        'led', 
        'readCount', 
        'rank',      // Added
        'rpu',       // Added
        'sequence'
    ];

    tabData.columnFilters = {
      sequenceID: '',
      sequence: '',
      cluster: '',
      readCount: '',
      rank: '',          // New
      rpu: '',           // New
      rankInCluster: '', // New
      led: ''            // New
    };
    
    tabData.tableData = data;
    tabData.filteredData = [...data];
    tabData.showResults = true;
    
    // Store the input filename for future use
    if (inputFile) {
      this.lastProcessedFile = inputFile;
    }
    
    // Ensure we're on the correct tab to see results
    if (this.selectedTabIndex() !== tabIndex) {
      this.selectedTabIndex.set(tabIndex);
    }
    
    this.updatePagination();
    
    console.log('Cluster results ready:', data.length, 'sequences, showing on tab', tabIndex);
    console.log('Input file stored:', this.lastProcessedFile);
  }

  private async loadPlotly(): Promise<void> {
    try {
      this.plotly = await import('plotly.js-dist-min');
    } catch (error) {
      console.error('Failed to load Plotly:', error);
    }
  }

  private async waitForPlotly(): Promise<void> {
    // If not on browser, return early
    if (!isPlatformBrowser(this.platformId)) return;
    
    // Wait for Plotly to load
    let attempts = 0;
    while (!this.plotly && attempts < 50) {
      await new Promise(resolve => setTimeout(resolve, 100));
      attempts++;
    }
    
    if (!this.plotly) {
      console.error('Plotly failed to load after waiting');
    }
  }

  async openMetaplotModal(event: { data: any, params: any }): Promise<void> {
    this.metaplotData = event.data;
    this.metaplotParams = event.params;
    this.showMetaplotModalFlag = true;
    
    // Wait for Plotly to load and modal to render
    await this.waitForPlotly();
    setTimeout(() => {
      this.createMetaplot();
    }, 100);
  }

  closeMetaplotModal(): void {
    this.showMetaplotModalFlag = false;
    this.metaplotData = null;
    this.metaplotParams = null;
  }

  async openKmerPlotModal(event: { data: any, params: any }): Promise<void> {
    this.kmerPlotData = event.data;
    this.kmerPlotParams = event.params;
    this.showKmerPlotModalFlag = true;
    
    // Wait for Plotly to load and modal to render
    await this.waitForPlotly();
    setTimeout(() => {
      this.createKmerPlot();
    }, 100);
  }

  closeKmerPlotModal(): void {
    this.showKmerPlotModalFlag = false;
    this.kmerPlotData = null;
    this.kmerPlotParams = null;
  }

  private createMetaplot(): void {
    if (!this.metaplotData || !this.metaplotParams || !this.plotly) {
      console.warn('Cannot create metaplot: Missing data, params, or Plotly not loaded');
      return;
    }

    // Prepare data for three separate line plots (matching R implementation)
    // R code: fa_clusterDiversity_metaplot creates 3 separate plots using plotly::subplot
    const clusterNumbers = this.metaplotData.map((d: any) => d.cluster);
    const totalSequences = this.metaplotData.map((d: any) => d.totalSequences);
    const totalReads = this.metaplotData.map((d: any) => d.totalReads);
    const averageLED = this.metaplotData.map((d: any) => d.averageLED);

    // Create 3 separate traces for subplot layout (matching R plotly::subplot approach)
    const trace1 = {
      x: clusterNumbers,
      y: totalSequences,
      type: 'scatter',
      mode: 'lines',
      name: this.metaplotParams.yaxes[0],
      line: { color: this.metaplotParams.line_colours[0], width: 2 },
      xaxis: 'x1',
      yaxis: 'y1'
    };

    const trace2 = {
      x: clusterNumbers,
      y: totalReads,
      type: 'scatter',
      mode: 'lines',
      name: this.metaplotParams.yaxes[1],
      line: { color: this.metaplotParams.line_colours[1], width: 2 },
      xaxis: 'x2',
      yaxis: 'y2'
    };

    const trace3 = {
      x: clusterNumbers,
      y: averageLED,
      type: 'scatter',
      mode: 'lines',
      name: this.metaplotParams.yaxes[2],
      line: { color: this.metaplotParams.line_colours[2], width: 2 },
      xaxis: 'x3',
      yaxis: 'y3'
    };

    // Layout mimicking R's plotly::subplot with 3 rows, 1 column, shareX = TRUE, titleY = TRUE
    const layout = {
      title: {
        text: `<b>${this.metaplotParams.plot_title}</b>`,
        font: { size: 16, family: 'Arial, sans-serif' },
        x: 0.5
      },
      // First subplot (top) - Total Sequences
      xaxis: { 
        domain: [0, 1],
        anchor: 'y1',
        showticklabels: false, // Hide x-axis labels for top subplot
        tickfont: { size: 12 }
      },
      yaxis: {
        title: {
          text: this.metaplotParams.yaxes[0],
          font: { 
            color: this.metaplotParams.line_colours[0],
            size: 14,
            family: 'Arial, sans-serif'
          }
        },
        domain: [0.7, 1],
        anchor: 'x1',
        tickfont: { 
          color: this.metaplotParams.line_colours[0],
          size: 12
        }
      },
      // Second subplot (middle) - Total Reads
      xaxis2: { 
        domain: [0, 1],
        anchor: 'y2',
        showticklabels: false, // Hide x-axis labels for middle subplot
        tickfont: { size: 12 }
      },
      yaxis2: {
        title: {
          text: this.metaplotParams.yaxes[1],
          font: { 
            color: this.metaplotParams.line_colours[1],
            size: 14,
            family: 'Arial, sans-serif'
          }
        },
        domain: [0.35, 0.65],
        anchor: 'x2',
        tickfont: { 
          color: this.metaplotParams.line_colours[1],
          size: 12
        }
      },
      // Third subplot (bottom) - Average LED
      xaxis3: { 
        title: {
          text: this.metaplotParams.xaxis,
          font: { size: 14, family: 'Arial, sans-serif' }
        },
        domain: [0, 1],
        anchor: 'y3',
        tickfont: { size: 12 }
      },
      yaxis3: {
        title: {
          text: this.metaplotParams.yaxes[2],
          font: { 
            color: this.metaplotParams.line_colours[2],
            size: 14,
            family: 'Arial, sans-serif'
          }
        },
        domain: [0, 0.3],
        anchor: 'x3',
        tickfont: { 
          color: this.metaplotParams.line_colours[2],
          size: 12
        }
      },
      height: 750,
      width: 950,
      margin: {
        l: 80,
        r: 80,
        t: 80,
        b: 80
      },
      showlegend: false, // R version doesn't show legend for line colors since axes are colored
      plot_bgcolor: 'white',
      paper_bgcolor: 'white'
    };

    const metaplotElement = document.getElementById('metaplot-container');
    if (metaplotElement) {
      // Create plot with config matching R's toImageButtonOptions
      const config = {
        toImageButtonOptions: {
          format: 'svg',
          filename: 'cluster_metaplots'
        },
        displayModeBar: true
      };
      
      this.plotly.newPlot(metaplotElement, [trace1, trace2, trace3], layout, config);
    }
  }

  private createKmerPlot(): void {
    if (!this.kmerPlotData || !this.kmerPlotParams || !this.plotly) {
      console.warn('Cannot create k-mer plot: Missing data, params, or Plotly not loaded');
      return;
    }

    console.log('Creating k-mer plot with backend integration');
    
    const inputFile = this.getInputFileName();
    if (!inputFile) {
      console.error('No input file available for k-mer analysis');
      this.showKmerError('No input file available. Please run clustering first.');
      return;
    }

    // Call backend API for k-mer analysis
    this.performKmerAnalysis(
      this.kmerPlotParams.clustersToPlot,
      this.kmerPlotParams.plotType.toLowerCase(),
      this.kmerPlotParams.kmerSize
    );
  }

  /**
   * Perform k-mer analysis using backend API with proper type safety
   */
  private performKmerAnalysis(clusters: number[], method: string, kSize: number): void {
    const requestBody: KmerAnalysisRequest = {
      input_path: this.getInputFileName(),
      selected_clusters: clusters,
      k_size: kSize,
      method: method as 'pca' | 'umap'
    };

    // Show loading state
    this.showKmerLoading();

    this.apiService.post('/api/v1/cluster-kmer-analysis', requestBody).subscribe({
      next: (response: any) => {
        if (response.status === 'ok' && response.data) {
          this.renderKmerPlot(response.data, method, kSize);
        } else {
          console.error('Invalid response format:', response);
          this.showKmerError('Invalid response from k-mer analysis');
        }
      },
      error: (error: any) => {
        console.error('K-mer analysis failed:', error);
        this.showKmerError(`K-mer analysis failed: ${error.error?.detail || error.message || 'Unknown error'}`);
      }
    });
  }

  /**
   * Get the current input file name for k-mer analysis
   */
  private getInputFileName(): string {
    // Try to get the filename from the current tab's results or from a stored property
    return this.lastProcessedFile || '';
  }

  /**
   * Render k-mer plot with real analysis data using proper type safety
   */
  private renderKmerPlot(data: KmerAnalysisData, method: string, kSize: number): void {
    try {
      // Extract unique clusters and assign colors using color service
      const clusters = data.coordinates.map((coord: any) => Number(coord.cluster));
      const uniqueSet = new Set(clusters);
      const uniqueClusters = Array.from(uniqueSet).sort((a, b) => (a as number) - (b as number)) as number[];
      const colorMap = ClusterColorService.createColorMap(uniqueClusters, this.getSelectedColorPalette());

      // Group coordinates by cluster
      const clusterMap = new Map();
      data.coordinates.forEach((coord: any) => {
        if (!clusterMap.has(coord.cluster)) {
          clusterMap.set(coord.cluster, {
            x: [],
            y: [],
            sequences: [],
            color: colorMap[coord.cluster],
            cluster: coord.cluster
          });
        }
        const clusterData = clusterMap.get(coord.cluster);
        clusterData.x.push(coord.x);
        clusterData.y.push(coord.y);
        clusterData.sequences.push(coord.sequence);
      });

      // Create Plotly traces with proper typing
      const traces: PlotlyTrace[] = Array.from(clusterMap.values()).map(clusterData => ({
        x: clusterData.x,
        y: clusterData.y,
        mode: 'markers' as const,
        type: 'scatter' as const,
        name: `Cluster ${clusterData.cluster}`,
        marker: {
          color: clusterData.color,
          size: 8,
          opacity: 0.7
        },
        text: clusterData.sequences.map((seq: string) => 
          `Cluster: ${clusterData.cluster}<br>Sequence: ${this.truncateSequence(seq, 20)}...`
        ),
        textposition: 'top center' as const,
        hovertemplate: '<b>%{text}</b><br>X: %{x:.3f}<br>Y: %{y:.3f}<extra></extra>'
      }));

      // Create layout with proper axis labels and typing
      const layout: PlotlyLayout = {
        title: `K-mer Analysis (${data.method_info.method}) - K=${kSize}`,
        xaxis: { 
          title: data.axis_labels.x,
          showgrid: true 
        },
        yaxis: { 
          title: data.axis_labels.y,
          showgrid: true 
        },
        showlegend: true,
        hovermode: 'closest',
        plot_bgcolor: 'white',
        paper_bgcolor: 'white',
        height: 600,
        width: 800,
        margin: { l: 80, r: 120, t: 80, b: 80 }
      };

      const kmerplotElement = document.getElementById('kmerplot-container');
      if (kmerplotElement) {
        const config = {
          toImageButtonOptions: { format: 'svg', filename: 'kmer_cluster' },
          displayModeBar: true
        };
        
        this.plotly.newPlot(kmerplotElement, traces, layout, config);
        console.log('K-mer plot rendered successfully with real data');
      }
    } catch (error) {
      console.error('Error rendering k-mer plot:', error);
      this.showKmerError('Failed to render k-mer plot');
    }
  }

  /**
   * Show loading state for k-mer plot
   */
  private showKmerLoading(): void {
    const kmerplotElement = document.getElementById('kmerplot-container');
    if (kmerplotElement && this.plotly) {
      const loadingLayout = {
        title: 'Loading K-mer Analysis...',
        xaxis: { title: 'Loading...', showgrid: false },
        yaxis: { title: 'Loading...', showgrid: false },
        showlegend: false,
        plot_bgcolor: 'white',
        paper_bgcolor: 'white',
        height: 600,
        width: 800,
        annotations: [{
          text: 'Performing k-mer analysis...',
          showarrow: false,
          xref: 'paper',
          yref: 'paper',
          x: 0.5,
          y: 0.5,
          xanchor: 'center',
          yanchor: 'middle',
          font: { size: 16, color: '#666' }
        }]
      };
      
      this.plotly.newPlot(kmerplotElement, [], loadingLayout, {});
    }
  }

  /**
   * Show error state for k-mer plot
   */
  private showKmerError(errorMessage: string): void {
    const kmerplotElement = document.getElementById('kmerplot-container');
    if (kmerplotElement && this.plotly) {
      const errorLayout = {
        title: 'K-mer Analysis Error',
        xaxis: { title: 'Error', showgrid: false },
        yaxis: { title: 'Error', showgrid: false },
        showlegend: false,
        plot_bgcolor: 'white',
        paper_bgcolor: 'white',
        height: 600,
        width: 800,
        annotations: [{
          text: errorMessage,
          showarrow: false,
          xref: 'paper',
          yref: 'paper',
          x: 0.5,
          y: 0.5,
          xanchor: 'center',
          yanchor: 'middle',
          font: { size: 14, color: '#d32f2f' }
        }]
      };
      
      this.plotly.newPlot(kmerplotElement, [], errorLayout, {});
    }
  }





  // Pagination methods (copied from start component)
  updatePagination(): void {
    const tabData = this.getCurrentTabData();
    tabData.totalPages = Math.ceil(tabData.filteredData.length / tabData.pageSize);
    if (tabData.currentPage > tabData.totalPages) {
      tabData.currentPage = Math.max(1, tabData.totalPages);
    }
    this.updatePaginatedData();
  }

  updatePaginatedData(): void {
    const tabData = this.getCurrentTabData();
    const startIndex = (tabData.currentPage - 1) * tabData.pageSize;
    const endIndex = startIndex + tabData.pageSize;
    tabData.paginatedData = tabData.filteredData.slice(startIndex, endIndex);
  }

  onPageSizeChange(): void {
    const tabData = this.getCurrentTabData();
    tabData.currentPage = 1;
    this.updatePagination();
  }

  previousPage(): void {
    const tabData = this.getCurrentTabData();
    if (tabData.currentPage > 1) {
      tabData.currentPage--;
      this.updatePaginatedData();
    }
  }

  nextPage(): void {
    const tabData = this.getCurrentTabData();
    if (tabData.currentPage < tabData.totalPages) {
      tabData.currentPage++;
      this.updatePaginatedData();
    }
  }

  goToPage(page: number): void {
    const tabData = this.getCurrentTabData();
    tabData.currentPage = page;
    this.updatePaginatedData();
  }

  getPageNumbers(): number[] {
    const tabData = this.getCurrentTabData();
    const pages: number[] = [];
    const maxVisible = 5;
    let startPage = Math.max(1, tabData.currentPage - Math.floor(maxVisible / 2));
    let endPage = Math.min(tabData.totalPages, startPage + maxVisible - 1);
    
    if (endPage - startPage < maxVisible - 1) {
      startPage = Math.max(1, endPage - maxVisible + 1);
    }
    
    for (let i = startPage; i <= endPage; i++) {
      pages.push(i);
    }
    return pages;
  }

  getShowingStart(): number {
    const tabData = this.getCurrentTabData();
    return tabData.filteredData.length === 0 ? 0 : (tabData.currentPage - 1) * tabData.pageSize + 1;
  }

  getShowingEnd(): number {
    const tabData = this.getCurrentTabData();
    return Math.min(tabData.currentPage * tabData.pageSize, tabData.filteredData.length);
  }

  onSearch(): void {
    this.applyFilters();
  }

  applyColumnFilter(): void {
    this.applyFilters();
  }

  applyFilters(): void {
    const tabData = this.getCurrentTabData();
    tabData.filteredData = tabData.tableData.filter(row => {
      // Global search
      const matchesSearch = !tabData.searchText || 
        Object.values(row).some(val => 
          String(val).toLowerCase().includes(tabData.searchText.toLowerCase())
        );
      
      // Dynamic column filters - check all columns that exist in columnFilters
      const matchesColumns = Object.keys(tabData.columnFilters).every(columnKey => {
        const filterValue = tabData.columnFilters[columnKey];
        if (!filterValue) return true; // No filter applied for this column
        
        const rowValue = row[columnKey];
        if (rowValue === undefined || rowValue === null) return false;
        
        // Handle different data types appropriately
        if (typeof rowValue === 'string') {
          return String(rowValue).toLowerCase().includes(filterValue.toLowerCase());
        } else {
          return String(rowValue).includes(filterValue);
        }
      });
      
      return matchesSearch && matchesColumns;
    });
    
    tabData.currentPage = 1;
    this.updatePagination();
  }

  /**
   * Get the selected color palette for cluster visualization
   */
  private getSelectedColorPalette(): string {
    // Get from k-mer plot parameters or use default
    return this.kmerPlotParams?.colour_palette_disc || 'Dark2';
  }

  /**
   * Truncate sequence for display purposes
   */
  private truncateSequence(sequence: string, maxLength: number): string {
    if (sequence.length <= maxLength) {
      return sequence;
    }
    return sequence.substring(0, maxLength);
  }

  /**
   * Enhanced color management using the color service
   */
  private getClusterColors(clusters: number[], palette: string): string[] {
    return ClusterColorService.assignClusterColors(clusters, palette).map(cc => cc.color);
  }

  /**
   * Get table CSS class based on current data type
   */
  getTableClass(): string {
    const currentTabIndex = this.selectedTabIndex();
    const tabData = this.tabResults.get(currentTabIndex);
    
    if (!tabData || !tabData.displayedColumns) {
      return '';
    }
    
    // Check if it's cluster results (has sequenceID column)
    if (tabData.displayedColumns.includes('sequenceID')) {
      return 'cluster-results';
    }
    
    // Check if it's diversity results (has totalSequences column)
    if (tabData.displayedColumns.includes('totalSequences')) {
      return 'diversity-results';
    }
    
    return '';
  }
}
