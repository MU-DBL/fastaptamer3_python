import { Component, Input, OnInit, OnChanges, SimpleChanges } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MatTableModule } from '@angular/material/table';

export interface ColumnDefinition {
  key: string;
  label: string;
  width?: string;
  exect_match?:boolean;
  filterable?: boolean;
}

export interface TableConfig {
  columns: ColumnDefinition[];
  initialPageSize?: number;
  pageSizeOptions?: number[];
}


@Component({
  selector: 'app-table',
  standalone: true,
  imports: [CommonModule, FormsModule, MatTableModule],
  templateUrl: './table.html',
  styleUrl: './table.scss',
})


export class Table implements OnInit, OnChanges {
  @Input() data: any[] = [];
  @Input() config: TableConfig = { columns: [] };
  
  // Table state
  filteredData: any[] = [];
  paginatedData: any[] = [];
  displayedColumns: string[] = [];
  
  // Pagination state
  currentPage: number = 1;
  pageSize: number = 10;
  totalPages: number = 0;
  pageSizeOptions: number[] = [10, 25, 50, 100];
  
  // Filter state
  searchText: string = '';
  columnFilters: { [key: string]: string } = {};

  ngOnInit(): void {
    this.initialize();
  }

  ngOnChanges(changes: SimpleChanges): void {
    if (changes['data'] || changes['config']) {
      this.initialize();
    }
  }

  private initialize(): void {
    // Set up columns
    this.displayedColumns = this.config.columns.map(col => col.key);
    
    // Initialize column filters
    this.columnFilters = {};
    this.config.columns.forEach(col => {
      if (col.filterable !== false) {
        this.columnFilters[col.key] = '';
      }
    });
    
    // Set page size options
    if (this.config.pageSizeOptions) {
      this.pageSizeOptions = this.config.pageSizeOptions;
    }
    
    // Set initial page size
    if (this.config.initialPageSize) {
      this.pageSize = this.config.initialPageSize;
    }
    
    // Initialize data
    this.filteredData = [...this.data];
    this.updatePagination();
  }

  // Pagination methods
  updatePagination(): void {
    this.totalPages = Math.ceil(this.filteredData.length / this.pageSize);
    if (this.currentPage > this.totalPages) {
      this.currentPage = Math.max(1, this.totalPages);
    }
    this.updatePaginatedData();
  }

  updatePaginatedData(): void {
    const startIndex = (this.currentPage - 1) * this.pageSize;
    const endIndex = startIndex + this.pageSize;
    this.paginatedData = this.filteredData.slice(startIndex, endIndex);
  }

  onPageSizeChange(): void {
    this.currentPage = 1;
    this.updatePagination();
  }

  previousPage(): void {
    if (this.currentPage > 1) {
      this.currentPage--;
      this.updatePaginatedData();
    }
  }

  nextPage(): void {
    if (this.currentPage < this.totalPages) {
      this.currentPage++;
      this.updatePaginatedData();
    }
  }

  goToPage(page: number): void {
    this.currentPage = page;
    this.updatePaginatedData();
  }

  getPageNumbers(): number[] {
    const pages: number[] = [];
    const maxVisible = 5;
    let startPage = Math.max(1, this.currentPage - Math.floor(maxVisible / 2));
    let endPage = Math.min(this.totalPages, startPage + maxVisible - 1);
    
    if (endPage - startPage < maxVisible - 1) {
      startPage = Math.max(1, endPage - maxVisible + 1);
    }
    
    for (let i = startPage; i <= endPage; i++) {
      pages.push(i);
    }
    return pages;
  }

  getShowingStart(): number {
    return this.filteredData.length === 0 ? 0 : (this.currentPage - 1) * this.pageSize + 1;
  }

  getShowingEnd(): number {
    return Math.min(this.currentPage * this.pageSize, this.filteredData.length);
  }

  // Filter methods
  onSearch(): void {
    this.applyFilters(false);
  }

  applyColumnFilter(isExactMatch: boolean = false): void {
    this.applyFilters(isExactMatch);
  }

  applyFilters(exect_match: boolean): void {
    this.filteredData = this.data.filter(row => {
      // Global search
      const matchesSearch = !this.searchText || 
        Object.values(row).some(val => 
          String(val).toLowerCase().includes(this.searchText.toLowerCase())
        );
      
      // Column filters
      const matchesColumns = Object.keys(this.columnFilters).every(key => {
        const filterValue = this.columnFilters[key];
        if (!filterValue) return true;

        const cellValue = String(row[key]).toLowerCase();
        const filterText = String(filterValue).toLowerCase();

        // --- LOGIC CHANGE HERE ---
        if (exect_match) {
           return cellValue === filterText; 
        } else {
           // Partial match (standard behavior)
           return cellValue.includes(filterText);
        }
      }); 
      return matchesSearch && matchesColumns;
    });
    
    this.currentPage = 1;
    this.updatePagination();
  }

  // Helper method to check if a column is filterable
  isColumnFilterable(columnKey: string): boolean {
    const column = this.config.columns.find(col => col.key === columnKey);
    return column?.filterable !== false;
  }

  // Helper method to get column label
  getColumnLabel(columnKey: string): string {
    const column = this.config.columns.find(col => col.key === columnKey);
    return column?.label || columnKey;
  }

  // Add these helper methods to your component

  formatCellValue(value: any, columnKey: string): string {
    if (value === null || value === undefined) {
      return '-';
    }
    
    // Special formatting for sequences
    if (columnKey === 'seqs' && typeof value === 'string' && value.length > 50) {
      return value; // Could add formatting here if needed
    }
    
    // Format numbers with commas
    if (typeof value === 'number' && columnKey !== 'rank') {
      return value.toLocaleString();
    }
    
    return String(value);
  }

  getCellTooltip(row: any, columnKey: string): string {
    const value = row[columnKey];
    
    if (value === null || value === undefined) {
      return '';
    }
    
    // Only show tooltip for long content
    if (typeof value === 'string' && value.length > 30) {
      return value;
    }
    
    return '';
  }

  isTruncated(value: any, maxWidth: string | undefined): boolean {
    if (!value || typeof value !== 'string') {
      return false;
    }
    
    // Simple heuristic: if value is long and width is constrained
    return value.length > 50 && !!maxWidth;
  }
}
