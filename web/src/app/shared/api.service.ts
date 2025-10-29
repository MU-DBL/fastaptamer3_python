import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { environment } from '../../environments/environment';

@Injectable({
  providedIn: 'root'
})
export class ApiService {
  private baseUrl = environment.apiUrl;

  constructor(private http: HttpClient) {}

  // Upload file
  uploadFile(file: File): Observable<any> {
    const formData = new FormData();
    formData.append('file', file);
    
    return this.http.post(`${this.baseUrl}/upload`, formData);
  }

  // Preprocess
  preprocess(params: {
    input_path: string;
    const5p: string;
    const3p: string;
    min_length: number;
    max_length: number;
    max_error: number;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/preprocess`, params);
  }

  // Count
  count(params: {
    input_path: string;
    reverseComplement: boolean;
    scaling_factor: number;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/count`, params);
  }

  // Download file
  downloadFile(filename: string): Observable<Blob> {
    return this.http.get(`${this.baseUrl}/download/${filename}`, {
      responseType: 'blob'
    });
  }

  // List files (optional, for future use)
  listFiles(): Observable<any> {
    return this.http.get(`${this.baseUrl}/files/`);
  }

  // Delete file (optional, for future use)
  deleteFile(filename: string): Observable<any> {
    return this.http.delete(`${this.baseUrl}/delete/${filename}`);
  }
}
