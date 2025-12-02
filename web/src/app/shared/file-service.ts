import { inject, Injectable } from '@angular/core';
import { ApiService } from './api.service';
import { Observable } from 'rxjs';


export enum ColumnName {
  ID = 'ID',
  RANK = 'Rank',
  READS = 'Reads',
  RPU = 'RPU',
  SEQUENCES = 'sequences',
  LENGTH = 'length',

  // Cluster columns
  CLUSTER = 'Cluster',
  RANK_IN_CLUSTER = 'RankInCluster',
  LED = 'LED',

  ORIGINAL_ID = 'OriginalID',

  // Diversity analysis columns
  TOTAL_SEQUENCES = 'TotalSequences',
  TOTAL_READS = 'TotalReads',
  TOTAL_RPU = 'TotalRPU',
  AVERAGE_LED = 'AverageLED',
  SID = 'SID'
}


@Injectable({
  providedIn: 'root'
})
export class FileService {

  private apiService = inject(ApiService);

  downloadFile(filename: string): void {
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
        alert(`Download failed: ${errorMsg}`);
      }
    });
  }

  parseClusterFile(blob: Blob, filename: string): Observable<any[]> {
    return new Observable(observer => {
      const reader = new FileReader();

      reader.onload = (e: any) => {
        const text = e.target.result;
        const parsedData = this.parseResultFile(text, filename);
        observer.next(parsedData);
        observer.complete();
      };

      reader.onerror = (error) => {
        observer.error(error);
      };
      reader.readAsText(blob);
    });
  }

  parseResultFile(content: string, filename: string): any[] {
    const data: any[] = [];

    if (filename.endsWith('.fasta') || filename.endsWith('.fa')) {
      return this.parseFastaFile(content);
    }else if( filename.endsWith('.csv') || filename.endsWith('.tsv') || filename.endsWith('.txt')) {
      return this.parseCsvFile(content);
    }
    console.warn('Unsupported file format for parsing:', filename);
    return data;
  }

  private parseFastaFile(content: string): any[] {

    const data: any[] = [];
    const fastaEntries = this.parseFasta(content);

    // 1. Create a Lookup Map for case-insensitive matching
    // Result: { 'rank': ColumnName.RANK, 'reads': ColumnName.READS, ... }
    const columnMap = new Map<string, ColumnName>();

    Object.values(ColumnName).forEach(val => {
      columnMap.set(val.toString().toLowerCase(), val);
    });

    fastaEntries.forEach((entry, index) => {
      // Initialize with defaults
      const row: any = {
        [ColumnName.ID]: entry.header,
        [ColumnName.SEQUENCES]: entry.sequence,
        [ColumnName.RANK]: index + 1,
        [ColumnName.READS]: 0,
        [ColumnName.RPU]: 0.0
      };

      const parts = entry.header.split(';');
      console.log('Parsing FASTA header parts:', parts);

      for (const part of parts) {
        if (!part.includes('=')) continue;

        const [rawKey, rawValue] = part.split('=', 2);
        const fileKey = rawKey.trim().toLowerCase();
        const valueString = rawValue.trim();

        // 3. Dynamic Lookup: Check if fileKey exists in our map
        const targetColumn = columnMap.get(fileKey);

        if (targetColumn) {
          // 4. Parse value based on target column type
          row[targetColumn] = this.parseValueByType(targetColumn, valueString);
        }
      }
      data.push(row);
    })
    console.log('Parsing FASTA  data:', data);
    return data;
  }

  private parseValueByType(column: ColumnName, value: string): any {
    const intColumns = [
      ColumnName.RANK,
      ColumnName.READS,
      ColumnName.CLUSTER,
      ColumnName.RANK_IN_CLUSTER,
      ColumnName.LED,
      ColumnName.TOTAL_SEQUENCES,
      ColumnName.TOTAL_READS,
      ColumnName.TOTAL_RPU,
      ColumnName.AVERAGE_LED,
      ColumnName.SID,
    ];

    // Define which columns are Floats
    const floatColumns = [
      ColumnName.RPU,
      ColumnName.TOTAL_RPU,
      ColumnName.AVERAGE_LED
    ];

    try {
      if (intColumns.includes(column)) {
        return parseInt(value, 10);
      }
      if (floatColumns.includes(column)) {
        return parseFloat(value);
      }
    } catch (e) {
      console.warn(`Error parsing number for ${column}`);
    }
    return value;
  }

  private parseFasta(content: string): { header: string; sequence: string }[] {
    const entries: { header: string; sequence: string }[] = [];
    const lines = content.split('\n');
    let currentEntry: { header: string; sequence: string } | null = null;

    for (const line of lines) {
      const trimmedLine = line.trim();
      if (trimmedLine.startsWith('>')) {
        if (currentEntry) {
          entries.push(currentEntry);
        }
        currentEntry = {
          header: trimmedLine.substring(1),
          sequence: ''
        };
      } else if (currentEntry && trimmedLine) {
        currentEntry.sequence += trimmedLine;
      }
    }

    if (currentEntry) {
      entries.push(currentEntry);
    }

    return entries;
  }

  private parseCsvFile(content: string): any[] {
    const data: any[] = [];

    // 1. Split content into lines and handle different line endings
    const lines = content.split(/\r?\n/).filter(line => line.trim().length > 0);

    if (lines.length === 0) {
      console.warn('CSV file is empty');
      return data;
    }

    // 2. Parse header row
    const headerLine = lines[0];
    const headers = this.parseCsvLine(headerLine);

    // 3. Create column mapping (case-insensitive)
    const columnMap = new Map<string, ColumnName>();
    Object.values(ColumnName).forEach(val => {
      columnMap.set(val.toString().toLowerCase(), val);
    });

    // 4. Map CSV headers to ColumnName enum
    const headerMapping: (ColumnName | null)[] = headers.map(header => {
      const normalizedHeader = header.trim().toLowerCase();
      return columnMap.get(normalizedHeader) || null;
    });

    // 5. Parse data rows
    for (let i = 1; i < lines.length; i++) {
      const values = this.parseCsvLine(lines[i]);

      // Skip rows that don't have enough columns
      if (values.length === 0) continue;

      const row: any = {};

      // 6. Map each value to its corresponding column
      values.forEach((value, colIndex) => {
        const targetColumn = headerMapping[colIndex];

        if (targetColumn) {
          row[targetColumn] = this.parseValueByType(targetColumn, value);
        }
      });

      // 7. Add default values for missing required columns
      if (!row[ColumnName.RANK]) {
        row[ColumnName.RANK] = i; // Use line number as rank if not provided
      }
      if (!row[ColumnName.READS]) {
        row[ColumnName.READS] = 0;
      }
      if (!row[ColumnName.RPU]) {
        row[ColumnName.RPU] = 0.0;
      }

      data.push(row);
    }

    console.log('Parsed CSV data:', data);
    return data;
  }

  /**
   * Parses a single CSV line, handling quoted fields and escaped quotes
   * Example: 'field1,"field2,with,commas","field3""with""quotes"'
   */
  private parseCsvLine(line: string): string[] {
    const result: string[] = [];
    let currentField = '';
    let insideQuotes = false;

    for (let i = 0; i < line.length; i++) {
      const char = line[i];
      const nextChar = line[i + 1];

      if (char === '"') {
        if (insideQuotes && nextChar === '"') {
          // Escaped quote: "" becomes "
          currentField += '"';
          i++; // Skip next quote
        } else {
          // Toggle quote state
          insideQuotes = !insideQuotes;
        }
      } else if (char === ',' && !insideQuotes) {
        // Field separator
        result.push(currentField.trim());
        currentField = '';
      } else {
        currentField += char;
      }
    }

    // Add the last field
    result.push(currentField.trim());
    return result;
  }
}



