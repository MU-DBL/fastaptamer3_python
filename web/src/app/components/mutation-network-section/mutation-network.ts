import {
  Component,
  inject,
  signal,
  ViewChild,
  ElementRef,
  AfterViewInit,
  OnDestroy,
  ChangeDetectorRef,
} from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { FileUploadResult, Upload } from '../../common/upload/upload';
import { ApiService } from '../../../shared/api.service';
import {
  MutationNetworkService,
  MutationNetworkGraphData,
} from './mutation-network.service';
import { switchMap, tap, catchError, finalize, timeout } from 'rxjs/operators';
import { of } from 'rxjs';
import { TimeoutError } from 'rxjs';
import * as d3 from 'd3';

/** Node with D3 simulation coordinates (mutated by forceSimulation). */
interface SimNode {
  id: string;
  sequence: string;
  index: number;
  x?: number;
  y?: number;
  fx?: number | null;
  fy?: number | null;
}

/** Link with source/target as node refs (used by D3). */
interface SimLink {
  source: SimNode;
  target: SimNode;
  cost: number;
}

@Component({
  selector: 'app-mutation-network',
  imports: [
    CommonModule,
    FormsModule,
    Upload,
    ...MATERIAL_IMPORTS,
  ],
  templateUrl: './mutation-network.html',
  styleUrl: './mutation-network.scss',
  standalone: true,
})
export class MutationNetwork implements AfterViewInit, OnDestroy {
  @ViewChild('graphContainer', { static: false }) graphContainer!: ElementRef<HTMLDivElement>;

  private apiService = inject(ApiService);
  private mutationService = inject(MutationNetworkService);
  private cdr = inject(ChangeDetectorRef);

  /** Reused from Count: same upload API (ApiService.uploadFile). */
  selectedFile: File | null = null;
  savedFileName = '';
  fileName = 'FASTA file';
  uploadComplete = false;

  startSequence = '';
  endSequence = '';
  maxDistance = 1;
  outputFormat: 'csv' | 'tsv' = 'csv';

  isProcessing = signal(false);
  processedFileName = signal('');
  graphData = signal<MutationNetworkGraphData | null>(null);
  errorMessage = signal<string | null>(null);

  /** Column IDs for path table (matches R output: From_Sequence, To_Sequence, Transition_Cost). */
  pathTableColumns = ['fromSequence', 'toSequence', 'transitionCost'] as const;

  private svg: d3.Selection<SVGSVGElement, unknown, null, undefined> | null = null;
  private resizeObserver: ResizeObserver | null = null;

  /** Start enabled only when file uploaded and start/end sequences provided (like Count: disabled until ready). */
  get canStart(): boolean {
    return (
      this.uploadComplete &&
      !!this.savedFileName?.trim() &&
      !!this.startSequence?.trim() &&
      !!this.endSequence?.trim() &&
      !this.isProcessing()
    );
  }

  onFileSelected(_result: FileUploadResult): void {
    this.selectedFile = _result.file;
    this.processedFileName.set('');
    this.graphData.set(null);
    this.errorMessage.set(null);
  }

  onUploadComplete(result: FileUploadResult): void {
    if (result.uploadComplete && result.savedFileName) {
      this.uploadComplete = true;
      this.savedFileName = result.savedFileName;
      this.errorMessage.set(null);
    } else if (result.error) {
      this.uploadComplete = false;
      this.errorMessage.set(result.error);
    }
  }

  onStart(): void {
    if (!this.canStart) return;

    this.isProcessing.set(true);
    this.processedFileName.set('');
    this.graphData.set(null);
    this.errorMessage.set(null);

    const params = {
      input_path: this.savedFileName,
      start_node: this.startSequence.trim(),
      end_node: this.endSequence.trim(),
      max_cost: this.maxDistance,
      output_format: this.outputFormat,
    };

    const REQUEST_TIMEOUT_MS = 120_000; // 2 minutes – backend can be slow for large files

    this.mutationService
      .run(params)
      .pipe(
        timeout(REQUEST_TIMEOUT_MS),
        tap((response) => {
          if (response.status === 'ok' && response.result) {
            this.processedFileName.set(response.result);
          }
        }),
        switchMap((response) => {
          if (response.status === 'ok' && response.result) {
            return this.apiService.downloadFile(response.result).pipe(
              timeout(30_000),
              tap((blob) => this.handleDownloadedResult(blob, response.result))
            );
          }
          return of(null);
        }),
        catchError((error) => {
          if (error instanceof TimeoutError) {
            this.errorMessage.set(
              'Request timed out. The file may be too large or the server is busy. Try a smaller file or try again.'
            );
          } else {
            const msg = error.error?.detail || error.message || 'Mutation network failed';
            this.errorMessage.set(Array.isArray(msg) ? msg.join(' ') : msg);
          }
          return of(null);
        }),
        finalize(() => {
          this.isProcessing.set(false);
          this.cdr.markForCheck();
        })
      )
      .subscribe();
  }

  private handleDownloadedResult(blob: Blob, _filename: string): void {
    const reader = new FileReader();
    reader.onload = (e: ProgressEvent<FileReader>) => {
      const text = (e.target?.result as string) || '';
      const data = this.mutationService.parseResultCsv(text);
      this.graphData.set(data);
      this.cdr.markForCheck();
      // Delay so view updates and #graphContainer (inside @if) is in DOM before renderGraph
      setTimeout(() => this.renderGraph(), 50);
    };
    reader.readAsText(blob);
  }

  onDownload(): void {
    const filename = this.processedFileName();
    if (!filename) return;

    this.apiService.downloadFile(filename).subscribe({
      next: (blob) => {
        const url = window.URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        link.click();
        window.URL.revokeObjectURL(url);
      },
      error: (error) => {
        this.errorMessage.set(error.error?.detail || 'Download failed');
      },
    });
  }

  ngAfterViewInit(): void {
    this.setupResizeObserver();
  }

  ngOnDestroy(): void {
    this.resizeObserver?.disconnect();
    this.resizeObserver = null;
    if (typeof window !== 'undefined') {
      window.removeEventListener('resize', this._onWindowResize);
    }
  }

  private setupResizeObserver(): void {
    const el = this.graphContainer?.nativeElement;
    if (!el) return;
    // Some server-side or older environments don't provide ResizeObserver
    // (e.g. SSR). Guard access and fall back to `window.resize`.
    const RO = (globalThis as any).ResizeObserver;
    if (typeof RO !== 'undefined') {
      const ro = new RO(() => {
        if (this.graphData()) this.renderGraph();
      });
      this.resizeObserver = ro;
      ro.observe(el);
    } else if (typeof window !== 'undefined') {
      // Fallback for environments without ResizeObserver
      window.addEventListener('resize', this._onWindowResize);
    }
  }

  // Window resize handler fallback
  private _onWindowResize = () => {
    if (this.graphData()) this.renderGraph();
  };

  private renderGraph(): void {
    const container = this.graphContainer?.nativeElement;
    const data = this.graphData();
    if (!container || !data || data.nodes.length === 0) {
      if (container) d3.select(container).selectAll('*').remove();
      return;
    }

    d3.select(container).selectAll('*').remove();
    const width = container.clientWidth;
    const height = container.clientHeight;
    if (width <= 0 || height <= 0) return;

    const svg = d3
      .select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .attr('viewBox', [0, 0, width, height]);

    const margin = { top: 40, right: 120, bottom: 40, left: 120 };
    const innerWidth = width - margin.left - margin.right;
    const innerHeight = height - margin.top - margin.bottom;

    const g = svg
      .append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    const nodes: SimNode[] = data.nodes.map((n) => ({ ...n }));
    const links = data.links.map((l) => ({ ...l }));

    const nodeById = new Map(nodes.map((n) => [n.id, n]));
    const linkObjects: SimLink[] = links.map((l) => ({
      source: nodeById.get(l.source) ?? (nodes[0] as SimNode),
      target: nodeById.get(l.target) ?? (nodes[0] as SimNode),
      cost: l.cost,
    }));

    const simulation = d3
      .forceSimulation(nodes)
      .force(
        'link',
        d3.forceLink(linkObjects).id((d: d3.SimulationNodeDatum) => (d as SimNode).id).distance(80)
      )
      .force('charge', d3.forceManyBody().strength(-200))
      .force('center', d3.forceCenter(innerWidth / 2, innerHeight / 2))
      .force('x', d3.forceX(innerWidth / 2).strength(0.05))
      .force('y', d3.forceY(innerHeight / 2).strength(0.05));

    const link = g
      .append('g')
      .selectAll<SVGLineElement, SimLink>('line')
      .data(linkObjects)
      .join('line')
      .attr('stroke', '#999')
      .attr('stroke-opacity', 0.6)
      .attr('stroke-width', 2);

    const linkLabels = g
      .append('g')
      .selectAll<SVGTextElement, SimLink>('text')
      .data(linkObjects)
      .join('text')
      .attr('font-size', 10)
      .attr('fill', '#333')
      .text((d) => String(d.cost));

    const node = g
      .append('g')
      .selectAll<SVGCircleElement, SimNode>('circle')
      .data(nodes)
      .join('circle')
      .attr('r', 8)
      .attr('fill', '#1976d2')
      .attr('stroke', '#fff')
      .attr('stroke-width', 2);

    node.call(
      d3
        .drag<SVGCircleElement, SimNode>()
        .on('start', (event) => {
          if (!event.active) simulation.alphaTarget(0.3).restart();
          event.subject.fx = event.subject.x ?? 0;
          event.subject.fy = event.subject.y ?? 0;
        })
        .on('drag', (event) => {
          event.subject.fx = event.x;
          event.subject.fy = event.y;
        })
        .on('end', (event) => {
          if (!event.active) simulation.alphaTarget(0);
          event.subject.fx = null;
          event.subject.fy = null;
        })
    );

    const nodeLabels = g
      .append('g')
      .selectAll<SVGTextElement, SimNode>('text')
      .data(nodes)
      .join('text')
      .attr('font-size', 11)
      .attr('dx', 12)
      .attr('dy', 4)
      .attr('fill', '#333')
      .text((d) => (d.sequence.length > 20 ? d.sequence.slice(0, 20) + '…' : d.sequence));

    simulation.on('tick', () => {
      link
        .attr('x1', (d) => d.source.x ?? 0)
        .attr('y1', (d) => d.source.y ?? 0)
        .attr('x2', (d) => d.target.x ?? 0)
        .attr('y2', (d) => d.target.y ?? 0);
      linkLabels
        .attr('x', (d) => ((d.source.x ?? 0) + (d.target.x ?? 0)) / 2)
        .attr('y', (d) => ((d.source.y ?? 0) + (d.target.y ?? 0)) / 2);
      node.attr('cx', (d) => d.x ?? 0).attr('cy', (d) => d.y ?? 0);
      nodeLabels.attr('x', (d) => d.x ?? 0).attr('y', (d) => d.y ?? 0);
    });

    this.svg = svg;
  }
}
