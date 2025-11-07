/**
 * TypeScript interfaces for K-mer Analysis
 * 
 * Provides type safety for k-mer analysis data structures and API responses.
 */

export interface KmerCoordinate {
  x: number;
  y: number;
  cluster: number;
  sequence: string;
}

export interface ClusterSummary {
  cluster: number;
  count: number;
}

export interface KmerMethodInfo {
  method: 'PCA' | 'UMAP';
  explained_variance?: {
    PC1: number;
    PC2: number;
  };
}

export interface KmerAxisLabels {
  x: string;
  y: string;
}

export interface KmerParameters {
  k_size: number;
  method: string;
  total_sequences: number;
  total_clusters: number;
}

export interface KmerAnalysisData {
  coordinates: KmerCoordinate[];
  clusters: ClusterSummary[];
  method_info: KmerMethodInfo;
  axis_labels: KmerAxisLabels;
  parameters: KmerParameters;
}

export interface KmerAnalysisResponse {
  status: 'ok' | 'error';
  data?: KmerAnalysisData;
  error?: string;
}

export interface KmerAnalysisRequest {
  input_path: string;
  selected_clusters: number[];
  k_size: number;
  method: 'pca' | 'umap';
}

export interface PlotlyTrace {
  x: number[];
  y: number[];
  mode: 'markers';
  type: 'scatter';
  name: string;
  marker: {
    color: string;
    size: number;
    opacity: number;
  };
  text: string[];
  textposition: 'top center';
  hovertemplate: string;
}

export interface PlotlyLayout {
  title: string;
  xaxis: {
    title: string;
    showgrid: boolean;
  };
  yaxis: {
    title: string;
    showgrid: boolean;
  };
  showlegend: boolean;
  hovermode: string;
  plot_bgcolor: string;
  paper_bgcolor: string;
  height: number;
  width: number;
  margin: {
    l: number;
    r: number;
    t: number;
    b: number;
  };
}