import { Injectable, NgZone } from '@angular/core';
import { HttpClient } from '@angular/common/http';
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

  constructor(private http: HttpClient, private zone: NgZone) {}

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
        this.zone.run(() => {
          try {
            const data = JSON.parse(event.data);
            observer.next(data);
            if (data.stage === 'complete' || data.stage === 'error') {
              observer.complete();
              this.closeConnection(jobId);
            }
          } catch (error) {
            console.error('Error parsing SSE data:', error);
          }
        });
      };

      eventSource.onerror = (_error) => {
        if (eventSource.readyState === EventSource.CLOSED) {
          this.zone.run(() => {
            observer.error(new Error('SSE connection closed'));
            this.closeConnection(jobId);
          });
        }
        // readyState === CONNECTING means browser is auto-reconnecting — ignore
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

  cancelProcesses(): Observable<any> {
    return this.http.post(`${this.baseUrl}/cancel`, {});
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
    input_paths: string[];
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

  // Cluster
  clusterLed(params: {
    input_path: string;
    output_format: string;
    min_reads: number;
    max_led: number;
    total_clusters: number;
    keep_nc: boolean;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/clusterled`, params);
  }

  // Cluster Diversity
  clusterDiversity(params: {
    input_path: string;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/clusterdiversity`, params);
  }

  // K-mer Analysis
  clusterKmerAnalysis(params: {
    input_path: string;
    selected_clusters: number[];
    k_size: number;
    method: 'pca' | 'umap';
  }): Observable<{
    status: 'ok' | 'error';
    data?: any;
    error?: string;
  }> {
    return this.http.post<{
      status: 'ok' | 'error';
      data?: any;
      error?: string;
    }>(`${this.baseUrl}/cluster-kmer-analysis`, params);
  }

  // Generic post method for custom endpoints
  post(endpoint: string, body: any): Observable<any> {
    const url = endpoint.startsWith('/api/') ? `${this.baseUrl.replace('/api/v1', '')}${endpoint}` : `${this.baseUrl}${endpoint}`;
    return this.http.post(url, body);
  }


  // POST /clustermsa
  clusterMsa(params: {
    input_path: string;
    output_format: string;
    seq_type: string;
    cluster_selected: number;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/clustermsa`, params);
  }

  // POST /cluster-msa-entropy
  clusterMsaEntropy(params: {
    input_path: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/cluster-msa-entropy`, params);
  }

  // POST /cluster-msa-mutinfo
  clusterMsaMutInfo(params: {
    input_path: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/cluster-msa-mutinfo`, params);
  }

  // POST /cluster-phmm-simulate
  clusterPhmmSimulate(params: {
    input_path: string;
    num_sequences: number;
    sequence_length: number;
    output_format_phmm: string;
    output_format_simulation: string;
    pseudocount_method: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/cluster-phmm-simulate`, params);
  }

  getClusterList(params: { input_path: string }): Observable<any>  {
    return this.http.post(`${this.baseUrl}/cluster-list`, params);
  }

  getPositionEnrichment(params: {
    fadf_recluster_path: string;
    output_format: string;
    seq_type: string;
    cluster_selection?: number | null;
    cluster_column?: string | null;
  }): Observable<any>  {
    return this.http.post(`${this.baseUrl}/position-enrichment`, params);
  }

  // Recluster two populations
  recluster(params: {
    fadf1_cluster_path: string;
    fadf2_cluster_path: string;
    led_threshold: number;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/recluster`, params);
  }

  reclusterMulti(params: {
    fadf1_cluster_path: string;
    fadf2_cluster_path: string;
    fadf3_cluster_path: string;
    round1_label: string;
    round2_label: string;
    round3_label: string;
    led_threshold: number;
    output_format: string;
  }) {
    return this.http.post(`${this.baseUrl}/recluster-multi`, params);
  }

  // Get LED matrix between two populations
  getReclusterLedMatrix(params: {
    fadf1_cluster_path: string;
    fadf2_cluster_path: string;
    led_threshold?: number;
    use_parallel?: boolean;
    n_jobs?: number;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/recluster-led-matrix`, params);
  }

  // Motif Search
  motifSearch(params: {
    input_path: string;
    motif: string;
    highlight: boolean;
    partial: boolean;
    motif_type: string;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/motif-search`, params);
  }

  // Motif Omit
  motifOmit(params: {
    input_path: string;
    motif: string;
    partial: boolean;
    motif_type: string;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/motif-omit`, params);
  }

  // Motif Tracker
  motifTracker(params: {
    input_paths: string[];
    population_names: string[];
    query_list: string[];
    query_aliases?: string[];
    motif_type: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/motif-tracker`, params);
  }

  // Sequence Tracker
  sequenceTracker(params: {
    input_paths: string[];
    population_names: string[];
    query_list: string[];
    query_aliases?: string[];
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/sequence-tracker`, params);
  }

  // Motif Discovery
  motifDiscovery(params: {
    input_path: string;
    min_reads: number;
    length_range: number[];
    output_format: string;
    alphabet: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/motif-discovery`, params);
  }

  // Translate
  translate(params: {
    input_path: string;
    orf: number;
    converge: boolean;
    input_changes: Array<{Codon: string, Translation: string}> | null;
    translate_selection: string;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/translate`, params);
  }

  // Distance
  sequenceDistance(params: {
    input_path: string;
    query_sequence: string;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/sequence-distance`, params);
  }

  // Differential Analysis
  differentialAnalysis(params: {
    cond1_paths: string[];
    cond2_paths: string[];
    p_cutoff: number;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/differential-expression`, params);
  }

  // Mutation network
  mutationNetwork(params: {
    input_path: string;
    start_node: string;
    end_node: string;
    max_cost: number;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/mutation-network`, params);
  }

  // Data Merge
  dataMerge(params: {
    input_paths: string[];
    merge_type: string;
    output_format: string;
  }): Observable<any> {
    return this.http.post(`${this.baseUrl}/data-merge`, params);
  }

  // Sequence Persistence
  sequencePersistence(params: {
    merged_file_path: string;
  }): Observable<{
    status: string;
    data: Array<{ freq: number; seqCount: number }>;
  }> {
    return this.http.post<{
      status: string;
      data: Array<{ freq: number; seqCount: number }>;
    }>(`${this.baseUrl}/sequence-persistence`, params);
  }

  // UpSet Data
  upSetData(params: {
    merged_file_path: string;
    fasta_names?: string[];
  }): Observable<{
    status: string;
    sets: string[];
    set_sizes: { [key: string]: number };
    intersections: Array<{
      sets: string[];
      size: number;
      sequences: string[];
    }>;
    total_unique_sequences: number;
  }> {
    return this.http.post<{
      status: string;
      sets: string[];
      set_sizes: { [key: string]: number };
      intersections: Array<{
        sets: string[];
        size: number;
        sequences: string[];
      }>;
      total_unique_sequences: number;
    }>(`${this.baseUrl}/upset-data`, params);
  }
}



