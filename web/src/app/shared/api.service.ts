import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { environment } from '../../environments/environment';

export interface ProgressEvent {
  stage: string;
  message: string;
  progress: number | null;
  timestamp: number;
  data?: any;
}


@Injectable({
  providedIn: 'root'
})
export class ApiService {
  private baseUrl = environment.apiUrl;
  private eventSources = new Map<string, EventSource>();

  constructor(private http: HttpClient) {}

  startPreprocessing(params: {
    input_path: string;
    const5p: string;
    const3p: string;
    min_length: number;
    max_length: number;
    max_error: number;
    output_format: string;
  }): Observable<string> {
    return new Observable(observer => {
      this.http.post<{ result: string }>(`${this.baseUrl}/preprocess`, params)
        .subscribe({
          next: (response) => {
            observer.next(response.result);
            observer.complete();
          },
          error: (error) => observer.error(error)
        });
    });
  }

  subscribeToProgress(jobId: string): Observable<ProgressEvent> {
    return new Observable(observer => {
      const eventSource = new EventSource(`${this.baseUrl}/progress/${jobId}`);
      this.eventSources.set(jobId, eventSource);

      eventSource.onmessage = (event) => {
        try {
          const data = JSON.parse(event.data);
          observer.next(data);
          // console.log('✅ Parsed event:', data);
          // Complete stream if job is done
          if (data.stage === 'complete' || data.stage === 'error') {
            observer.complete();
            this.closeConnection(jobId);
          }
        } catch (error) {
          console.error('Error parsing SSE data:', error);
        }
      };

      eventSource.onerror = (error) => {
        console.error('SSE error:', error);
        observer.error(error);
        this.closeConnection(jobId);
      };

      // Cleanup function
      return () => {
        this.closeConnection(jobId);
      };
    });
  }

  private closeConnection(jobId: string): void {
    const eventSource = this.eventSources.get(jobId);
    if (eventSource) {
      eventSource.close();
      this.eventSources.delete(jobId);
    }
  }

  cancelJob(jobId: string): void {
    this.closeConnection(jobId);
  }

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

  // Recount
  recount(params: {
    input_path_1: string;
    input_path_2: string;
    scaling_factor: number;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/recount`, params);
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
