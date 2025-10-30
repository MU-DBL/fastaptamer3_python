import { Component, signal, WritableSignal } from '@angular/core';
import { Count } from '../count/count';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { Recount } from '../recount/recount';
import { Preprocess } from '../preprocess/preprocess';
import { ReadsPerRankChart } from '../count/charts/reads-per-rank-chart';
import { SeqLengthHistogram } from '../count/charts/seq-length-histogram';
import { AbundancePlot } from '../count/charts/abundance-plot';

@Component({
  selector: 'app-start',
  imports: [
    CommonModule,
    FormsModule,
    Count,
    Recount,
    Preprocess,
    ReadsPerRankChart,
    SeqLengthHistogram,
    AbundancePlot,
    ...MATERIAL_IMPORTS],
  templateUrl: './start.html',
  styleUrl: './start.scss'
})
export class Start {
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
    return this.tabResults.get(this.selectedTabIndex())?.displayedColumns ?? ['id', 'rank', 'reads', 'rpm', 'length', 'seqs'];
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
      id: '',
      rank: '',
      reads: '',
      rpm: '',
      length: '',
      seqs: ''
    };
  }

  set columnFilters(value: any) {
    const tabData = this.tabResults.get(this.selectedTabIndex());
    if (tabData) {
      tabData.columnFilters = value;
    }
  }

  // Modal state management
  showReadsPerRankModal = signal(false);
  showSeqLengthModal = signal(false);
  showAbundancePlotModal = signal(false);

  // Chart data
  readsPerRankData: any[] = [];
  seqLengthData: any = { unique: [], total: [] };
  abundancePlotData: any[] = [];

  // Chart customization parameters
  readsPerRankParams = {
    title: 'Read count for each rank',
    xAxisLabel: 'Ranks of unique sequences',
    yAxisLabel: 'Total reads per unique sequence',
    lineColor: '#87CEEB'
  };

  seqLengthParams = {
    title: 'Sequence-length histogram',
    xAxisLabel: 'Sequence length',
    yAxis1Label: 'Unique sequences',
    yAxis2Label: 'Read count',
    barOutline: '#000000',
    barFill: '#87CEEB',
    barFill2: '#FFA500'
  };

  abundancePlotParams = {
    title: 'Binned sequence abundance',
    xAxisLabel: 'No. Reads',
    yAxisLabel: 'Fraction of Population',
    barOutline: '#000000',
    barFill: '#87CEEB',
    colorLight: '#ADD8E6',
    colorDark: '#FF6B6B'
  };

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
        displayedColumns: ['id', 'rank', 'reads', 'rpm', 'length', 'seqs'],
        currentPage: 1,
        pageSize: 10,
        totalPages: 0,
        searchText: '',
        columnFilters: {
          id: '',
          rank: '',
          reads: '',
          rpm: '',
          length: '',
          seqs: ''
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

  onCountResultsReady(data: any[]): void {
    const tabIndex = this.selectedTabIndex(); // Count tab
    this.ensureTabData(tabIndex);
    const tabData = this.tabResults.get(tabIndex)!;
    
    tabData.tableData = data;
    tabData.filteredData = [...data];
    tabData.showResults = true;
    
    this.updatePagination();
  }

  onRecountResultsReady(data: any[]): void {
    const tabIndex = this.selectedTabIndex(); // Recount tab
    this.ensureTabData(tabIndex);
    const tabData = this.tabResults.get(tabIndex)!;
    
    tabData.tableData = data;
    tabData.filteredData = [...data];
    tabData.showResults = true;
    
    this.updatePagination();
  }

  // Modal event handlers from count component
  onShowReadsPerRankModal(event: { data: any[], params: any }): void {
    this.readsPerRankData = event.data;
    this.readsPerRankParams = { ...this.readsPerRankParams, ...event.params };
    this.showReadsPerRankModal.set(true);
  }

  onShowSeqLengthModal(event: { data: any, params: any }): void {
    this.seqLengthData = event.data;
    this.seqLengthParams = { ...this.seqLengthParams, ...event.params };
    this.showSeqLengthModal.set(true);
  }

  onShowAbundancePlotModal(event: { data: any[], params: any }): void {
    this.abundancePlotData = event.data;
    this.abundancePlotParams = { ...this.abundancePlotParams, ...event.params };
    this.showAbundancePlotModal.set(true);
  }

  closeReadsPerRankModal(): void {
    this.showReadsPerRankModal.set(false);
  }

  closeSeqLengthModal(): void {
    this.showSeqLengthModal.set(false);
  }

  closeAbundancePlotModal(): void {
    this.showAbundancePlotModal.set(false);
  }


  // Pagination methods
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
      
      // Column filters
      const matchesColumns = 
        (!tabData.columnFilters.id || String(row.id).toLowerCase().includes(tabData.columnFilters.id.toLowerCase())) &&
        (!tabData.columnFilters.rank || String(row.rank).includes(tabData.columnFilters.rank)) &&
        (!tabData.columnFilters.reads || String(row.reads).includes(tabData.columnFilters.reads)) &&
        (!tabData.columnFilters.rpm || String(row.rpm).includes(tabData.columnFilters.rpm)) &&
        (!tabData.columnFilters.length || String(row.length).includes(tabData.columnFilters.length)) &&
        (!tabData.columnFilters.seqs || String(row.seqs).toLowerCase().includes(tabData.columnFilters.seqs.toLowerCase()));
      
      return matchesSearch && matchesColumns;
    });
    
    tabData.currentPage = 1;
    this.updatePagination();
  }
}
