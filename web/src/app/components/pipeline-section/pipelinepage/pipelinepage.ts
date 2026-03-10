import { Component, inject, signal, OnDestroy, ChangeDetectorRef } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { lastValueFrom, map } from 'rxjs';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { Upload, FileUploadResult } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import { FileService } from '../../../shared/file-service';
import { GENETIC_CODES } from '../../../shared/constants';

interface PipelineStep {
  id: string;
  operationKey: string;
  params: Record<string, any>;  // keys are fixed; only values are editable
  status: 'idle' | 'running' | 'complete' | 'error';
  inputFile?: string;
  outputFile?: string;
  error?: string;
  progressLogs: string[];
  progressValue: number;
}

interface OperationDef {
  key: string;
  label: string;
  defaultParams: Record<string, any>;
  paramLabels: Record<string, string>;
  paramOptions?: Record<string, any[]>;  // fixed choices → renders as <mat-select>
  paramOptionLabels?: Record<string, Record<string, string>>;  // value → display label
  possibleNextSteps?: string[];
}

const OPERATION_DEFS: OperationDef[] = [
  {
    key: 'preprocess',
    label: 'Preprocess',
    defaultParams: {
      const5p: '',
      const3p: '',
      min_length: 10,
      max_length: 100,
      max_error: 0.005,
    },
    paramLabels: {
      const5p: "Constant 5' Region",
      const3p: "Constant 3' Region",
      min_length: 'Min Length',
      max_length: 'Max Length',
      max_error: 'Max Allowed Error (0-1)',
    },
    possibleNextSteps: ['count'],
  },
  {
    key: 'count',
    label: 'Count',
    defaultParams: {
      reverseComplement: false,
      scaling_factor: 1000000,
    },
    paramLabels: {
      reverseComplement: 'Reverse Complement',
      scaling_factor: 'Scaling Factor',
    },
    paramOptions: {
      reverseComplement: [false, true],
      scaling_factor: [1, 10, 100, 1000, 10000, 100000, 1000000],
    },
    possibleNextSteps: ['translate', 'motif-search', 'motif-omit', 'motif-discovery', 'mutation-network', 'sequence-distance', 'cluster'],
  },
  {
    key: 'translate',
    label: 'Translate',
    defaultParams: {
      orf: 1,
      converge: true,
      translate_selection: 'Standard',
    },
    paramLabels: {
      orf: 'Open Reading Frame',
      converge: 'Merge Identical Translations',
      translate_selection: 'Genetic Code',
    },
    paramOptions: {
      orf: [1, 2, 3],
      converge: [true, false],
      translate_selection: GENETIC_CODES,
    },
    possibleNextSteps: ['cluster','motif-discovery','mutation-network','sequence-distance'],
  },
  {
    key: 'motif-search',
    label: 'Motif Search',
    defaultParams: {
      motif: '',
      highlight: false,
      partial: false,
      motif_type: 'Nucleotide',
    },
    paramLabels: {
      motif: 'Motifs (comma-separated)',
      highlight: 'Highlight Matches in Sequence',
      partial: 'Allow Partial Matches (OR)',
      motif_type: 'Motif Type',
    },
    paramOptions: {
      highlight: [false, true],
      partial: [false, true],
      motif_type: ['Nucleotide', 'AminoAcid', 'String'],
    },
    possibleNextSteps: ['cluster','motif-discovery','mutation-network','sequence-distance'],
  },
  {
    key: 'motif-omit',
    label: 'Motif Omit',
    defaultParams: {
      motif: '',
      partial: false,
      motif_type: 'Nucleotide',
    },
    paramLabels: {
      motif: 'Motifs (comma-separated)',
      partial: 'Allow Partial Matches (OR)',
      motif_type: 'Motif Type',
    },
    paramOptions: {
      partial: [false, true],
      motif_type: ['Nucleotide', 'AminoAcid', 'String'],
    },
    possibleNextSteps: ['cluster','motif-discovery','mutation-network','sequence-distance'],
  },
  {
    key: 'motif-discovery',
    label: 'Motif Discovery',
    defaultParams: {
      min_reads: 10,
      length_min: 5,
      length_max: 10,
      alphabet: 'dna',
    },
    paramLabels: {
      min_reads: 'Min Reads',
      length_min: 'Min Motif Length',
      length_max: 'Max Motif Length',
      alphabet: 'Sequence Alphabet',
    },
    paramOptions: {
      alphabet: ['dna', 'protein'],
    },
    paramOptionLabels: {
      alphabet: { dna: 'Nucleotide', protein: 'AminoAcid' },
    },
  },
  {
    key: 'mutation-network',
    label: 'Mutation Network',
    defaultParams: {
      start_node: '',
      end_node: '',
      max_cost: 1,
    },
    paramLabels: {
      start_node: 'Start Sequence',
      end_node: 'End Sequence',
      max_cost: 'Max Edit Distance',
    }
  },
  {
    key: 'sequence-distance',
    label: 'Sequence Distance',
    defaultParams: {
      query_sequence: '',
    },
    paramLabels: {
      query_sequence: 'Query Sequence',
    }
  },
  {
    key: 'cluster',
    label: 'Cluster',
    defaultParams: {
      min_reads: 10,
      max_led: 7,
      total_clusters: 20,
      keep_nc: false,
    },
    paramLabels: {
      min_reads: 'Min Reads',
      max_led: 'Max LED',
      total_clusters: 'Max Clusters to Generate',
      keep_nc: 'Keep Non-Clustered',
    },
    paramOptions: {
      keep_nc: [false, true],
    },
    possibleNextSteps: ['cluster-msa','cluster-diversity'],
  },
  {
    key: 'cluster-msa',
    label: 'Cluster MSA',
    defaultParams: {
      cluster_selected: 1,
      seq_type: 'dna',
    },
    paramLabels: {
      cluster_selected: 'Cluster ID',
      seq_type: 'Sequence Type',
    },
    paramOptions: {
      seq_type: ['dna', 'protein'],
    },
    paramOptionLabels: {
      seq_type: { dna: 'Nucleotide', protein: 'AminoAcid' },
    },
    possibleNextSteps: ['cluster-phmm'],
  },
  {
    key: 'cluster-phmm',
    label: 'Cluster PHMM',
    defaultParams: {
      num_sequences: 100,
      sequence_length: 50,
    },
    paramLabels: {
      num_sequences: 'Number of Sequences to Simulate',
      sequence_length: 'Target Sequence Length',
    }
  },
  {
    key: 'cluster-diversity',
    label: 'Cluster Diversity',
    defaultParams: {},
    paramLabels: {},
  }
];

@Component({
  selector: 'app-pipelinepage',
  imports: [CommonModule, FormsModule, Upload, ...MATERIAL_IMPORTS],
  templateUrl: './pipelinepage.html',
  styleUrl: './pipelinepage.scss'
})

export class Pipelinepage implements OnDestroy {
  private apiService = inject(ApiService);
  private fileService = inject(FileService);
  private cdr = inject(ChangeDetectorRef);

  operationDefs = OPERATION_DEFS;
  steps: PipelineStep[] = [];
  initialFileName = '';
  isRunning = signal(false);
  private cancelled = false;

  onInitialUpload(result: FileUploadResult): void {
    if (this.initialFileName) {
      this.apiService.deleteFile(this.initialFileName).subscribe();
    }

    if (result.uploadComplete && result.savedFileName) {
      this.initialFileName = result.savedFileName;
    }
  }

  ngOnDestroy(): void {
    if (this.initialFileName) {
      this.apiService.deleteFile(this.initialFileName).subscribe();
    }
    
    this.steps.forEach(step => {
      if (step.inputFile) {
        this.apiService.deleteFile(step.inputFile).subscribe();
      }

      if (step.outputFile) {
        this.apiService.deleteFile(step.outputFile).subscribe();
      }
    });
  }

  // Preserve OPERATION_DEFS insertion order in keyvalue pipe
  readonly originalOrder = (): number => 0;

  addStep(): void {
    const defaultOp = OPERATION_DEFS[0];
    this.steps.push({
      id: Math.random().toString(36).substr(2, 9),
      operationKey: defaultOp.key,
      params: { ...defaultOp.defaultParams },
      status: 'idle',
      progressLogs: [],
      progressValue: 0
    });
  }

  removeStep(index: number): void {
    this.steps.splice(index, 1);
  }

  onOperationChange(step: PipelineStep): void {
    const def = OPERATION_DEFS.find(d => d.key === step.operationKey);
    if (def) {
      step.params = { ...def.defaultParams };
    }
    step.status = 'idle';
    step.outputFile = undefined;
    step.error = undefined;
    step.progressLogs = [];
    step.progressValue = 0;
  }

  getOperationLabel(key: string): string {
    return OPERATION_DEFS.find(d => d.key === key)?.label ?? key;
  }

  getParamLabel(opKey: string, paramKey: string): string {
    return OPERATION_DEFS.find(d => d.key === opKey)?.paramLabels?.[paramKey] ?? paramKey;
  }

  getParamOptions(opKey: string, paramKey: string): any[] | null {
    return OPERATION_DEFS.find(d => d.key === opKey)?.paramOptions?.[paramKey] ?? null;
  }

  getParamOptionLabel(opKey: string, paramKey: string, value: any): string {
    const label = OPERATION_DEFS.find(d => d.key === opKey)?.paramOptionLabels?.[paramKey]?.[value];
    return label ?? String(value);
  }

  getStepError(index: number): string | null {
    if (index === 0) return null;
    const prevDef = OPERATION_DEFS.find(d => d.key === this.steps[index - 1].operationKey);
    if (!prevDef?.possibleNextSteps) return null;
    const curKey = this.steps[index].operationKey;
    if (prevDef.possibleNextSteps.includes(curKey)) return null;
    const allowed = prevDef.possibleNextSteps.map(k => this.getOperationLabel(k)).join(', ');
    return `After "${prevDef.label}", next step must be one of: ${allowed}`;
  }

  canRun(): boolean {
    if (!this.initialFileName || this.steps.length === 0 || this.isRunning()) return false;
    return !this.steps.some((_, i) => this.getStepError(i) !== null);
  }

  cancelPipeline(): void {
    this.cancelled = true;
    this.apiService.cancelProcesses().subscribe();
    const runningStep = this.steps.find(s => s.status === 'running');
    if (runningStep) {
      runningStep.status = 'error';
      runningStep.error = 'Cancelled by user';
      this.cdr.detectChanges();
    }
    this.isRunning.set(false);
  }

  async runPipeline(): Promise<void> {
    this.cancelled = false;
    for (const step of this.steps) {
      step.status = 'idle';
      step.outputFile = undefined;
      step.error = undefined;
      step.progressLogs = [];
      step.progressValue = 0;
    }

    this.isRunning.set(true);
    let currentFile = this.initialFileName;

    for (const step of this.steps) {
      if (this.cancelled) break;
      step.status = 'running';
      step.inputFile = currentFile;
      this.cdr.detectChanges();
      try {
        const outputFile = await this.executeStep(step, currentFile);
        if (this.cancelled) break;
        if (!outputFile) throw new Error('Step returned no output file path');
        currentFile = outputFile;
        step.outputFile = currentFile;
        step.status = 'complete';
        this.cdr.detectChanges();
      } catch (err: any) {
        if (!this.cancelled) {
          step.status = 'error';
          step.error = err?.error?.detail || err?.message || 'Processing failed';
          this.cdr.detectChanges();
        }
        break;
      }
    }

    this.isRunning.set(false);
  }

  private async executeStep(step: PipelineStep, inputFile: string): Promise<string> {
    const p = (output_format = 'fasta', extra?: Record<string, any>) =>
      ({ input_path: inputFile, ...step.params, output_format, ...extra } as any);

    switch (step.operationKey) {
      case 'preprocess':
        return this.executePreprocess(step, inputFile);
      case 'count':
        return lastValueFrom(this.apiService.count(p()).pipe(map((r: any) => r.result)));
      case 'translate':
        return lastValueFrom(this.apiService.translate(p('fasta', { input_changes: null })).pipe(map((r: any) => r.result)));
      case 'motif-search':
        return lastValueFrom(this.apiService.motifSearch(p()).pipe(map((r: any) => r.result)));
      case 'motif-omit':
        return lastValueFrom(this.apiService.motifOmit(p()).pipe(map((r: any) => r.result)));
      case 'motif-discovery':
        return lastValueFrom(this.apiService.motifDiscovery({
          input_path: inputFile,
          min_reads: step.params['min_reads'],
          length_range: [step.params['length_min'], step.params['length_max']],
          output_format: 'csv',
          alphabet: step.params['alphabet'],
        }).pipe(map((r: any) => r.result)));
      case 'mutation-network':
        return lastValueFrom(this.apiService.mutationNetwork(p('csv')).pipe(map((r: any) => r.result)));
      case 'sequence-distance':
        return lastValueFrom(this.apiService.sequenceDistance(p('csv')).pipe(map((r: any) => r.result)));
      case 'cluster':
        return lastValueFrom(this.apiService.clusterLed(p()).pipe(map((r: any) => r.result)));
      case 'cluster-msa':
        return lastValueFrom(this.apiService.clusterMsa(p()).pipe(map((r: any) => r.result)));
      case 'cluster-phmm':
        return lastValueFrom(this.apiService.clusterPhmmSimulate(
          p('fasta', { output_format_phmm: 'txt', output_format_simulation: 'fasta' })
        ).pipe(map((r: any) => r.output_file_simulation)));
      case 'cluster-diversity':
        return lastValueFrom(this.apiService.clusterDiversity(p('csv')).pipe(map((r: any) => r.result)));
      default:
        throw new Error(`Unknown operation: ${step.operationKey}`);
    }
  }

  private executePreprocess(step: PipelineStep, inputFile: string): Promise<string> {
    return new Promise((resolve, reject) => {
      this.apiService.startPreprocessing({ input_path: inputFile, ...step.params, output_format: 'fasta' } as any)
        .subscribe({
          next: jobId => {
            this.apiService.subscribeToProgress(jobId).subscribe({
              next: event => {
                step.progressLogs.push(`[${event.stage}] ${event.message}`);
                if (event.progress != null) step.progressValue = event.progress;
                if (event.stage === 'complete' && event.data?.output_path) {
                  resolve(event.data.output_path);
                }
                if (event.stage === 'error') {
                  reject(new Error(event.message));
                }
              },
              error: reject
            });
          },
          error: reject
        });
    });
  }

  downloadResult(step: PipelineStep): void {
    if (!step.outputFile) return;
    this.fileService.downloadFile(step.outputFile);
  }
}
