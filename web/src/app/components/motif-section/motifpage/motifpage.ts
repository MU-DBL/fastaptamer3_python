import { Component, signal, WritableSignal } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MotifSearch } from '../motif-search/motif-search';
import { MotifOmit } from '../motif-omit/motif-omit';
import { MotifTracker } from '../motif-tracker/motif-tracker';
import { MotifDiscovery } from '../motif-discovery/motif-discovery';

@Component({
  selector: 'app-motifpage',
  imports: [
    CommonModule,
    FormsModule,
    MotifSearch,
    MotifOmit,
    MotifTracker,
    MotifDiscovery,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './motifpage.html',
  styleUrl: './motifpage.scss',
})
export class Motifpage {
  
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
           ['id', 'rank', 'reads', 'rpm', 'length', 'seqs'];
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
    return this.tabResults.get(this.selectedTabIndex())?.columnFilters ?? {};
  }

  set columnFilters(value: any) {
    const tabData = this.tabResults.get(this.selectedTabIndex());
    if (tabData) {
      tabData.columnFilters = value;
    }
  }

  // Handler for when results are ready from child components
  onResultsReady(data: any[]): void {
    const currentIndex = this.selectedTabIndex();
    
    // Initialize or update the current tab's data
    if (!this.tabResults.has(currentIndex)) {
      this.tabResults.set(currentIndex, {
        showResults: true,
        tableData: data,
        filteredData: data,
        paginatedData: data.slice(0, 10),
        displayedColumns: ['id', 'rank', 'reads', 'rpm', 'length', 'seqs'],
        currentPage: 1,
        pageSize: 10,
        totalPages: Math.ceil(data.length / 10),
        searchText: '',
        columnFilters: {}
      });
    } else {
      const tabData = this.tabResults.get(currentIndex)!;
      tabData.showResults = true;
      tabData.tableData = data;
      tabData.filteredData = data;
      tabData.paginatedData = data.slice(0, tabData.pageSize);
      tabData.totalPages = Math.ceil(data.length / tabData.pageSize);
      tabData.currentPage = 1;
    }
  }
}
