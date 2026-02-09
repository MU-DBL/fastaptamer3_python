import { Injectable } from '@angular/core';
import { ApiService } from '../../../shared/api.service';
import { Observable } from 'rxjs';

/** One step in the mutation path (backend CSV row). */
export interface MutationPathStep {
  fromSequence: string;
  toSequence: string;
  transitionCost: number;
}

/** Graph node for visualization. */
export interface MutationNetworkNode {
  id: string;
  sequence: string;
  index: number;
}

/** Graph link for visualization. */
export interface MutationNetworkLink {
  source: string;
  target: string;
  cost: number;
}

/** Parsed path and graph data from backend result. */
export interface MutationNetworkGraphData {
  steps: MutationPathStep[];
  nodes: MutationNetworkNode[];
  links: MutationNetworkLink[];
}

export interface MutationNetworkParams {
  input_path: string;
  start_node: string;
  end_node: string;
  max_cost: number;
  output_format: string;
}

/**
 * Service for Mutation Network feature.
 * Reuses ApiService.uploadFile and ApiService.downloadFile (same as Count).
 * Backend: POST /api/v1/mutation-network → returns result filename; download CSV to get path.
 */
@Injectable({
  providedIn: 'root'
})
export class MutationNetworkService {

  constructor(private api: ApiService) {}

  /** Run mutation network (same async pattern as Count: POST then optional download). */
  run(params: MutationNetworkParams): Observable<{ status: string; result: string }> {
    return this.api.mutationNetwork(params);
  }

  /** Download result file (reuses Count download API). */
  downloadFile(filename: string): Observable<Blob> {
    return this.api.downloadFile(filename);
  }

  /**
   * Parse CSV content from backend result into path steps and graph data.
   * Backend CSV columns: From_Sequence, To_Sequence, Transition_Cost.
   */
  parseResultCsv(content: string): MutationNetworkGraphData {
    const steps: MutationPathStep[] = [];
    const lines = content.split('\n').filter(line => line.trim());
    if (lines.length < 2) {
      return { steps: [], nodes: [], links: [] };
    }
    const sep = lines[0].includes('\t') ? '\t' : ',';
    const headers = lines[0].split(sep).map(h => h.replace(/^"|"$/g, '').trim().toLowerCase());
    const fromIdx = headers.findIndex(h => h === 'from_sequence') >= 0 ? headers.findIndex(h => h === 'from_sequence') : 0;
    const toIdx = headers.findIndex(h => h === 'to_sequence') >= 0 ? headers.findIndex(h => h === 'to_sequence') : 1;
    const costIdx = headers.findIndex(h => h === 'transition_cost') >= 0 ? headers.findIndex(h => h === 'transition_cost') : 2;
    const parseRow = (line: string): string[] => {
      const result: string[] = [];
      let inQuotes = false;
      let current = '';
      for (let i = 0; i < line.length; i++) {
        const c = line[i];
        if (c === '"') {
          inQuotes = !inQuotes;
        } else if (c === sep && !inQuotes) {
          result.push(current.trim());
          current = '';
        } else {
          current += c;
        }
      }
      result.push(current.trim());
      return result;
    };
    for (let i = 1; i < lines.length; i++) {
      const cols = parseRow(lines[i]);
      const fromSeq = cols[fromIdx] ?? '';
      const toSeq = cols[toIdx] ?? '';
      const cost = parseInt(cols[costIdx] ?? '0', 10) || 0;
      if (fromSeq && toSeq) {
        steps.push({ fromSequence: fromSeq, toSequence: toSeq, transitionCost: cost });
      }
    }
    const nodes: MutationNetworkNode[] = [];
    const links: MutationNetworkLink[] = [];
    const seen = new Set<string>();
    let idx = 0;
    if (steps.length > 0) {
      const first = steps[0];
      if (!seen.has(first.fromSequence)) {
        nodes.push({ id: first.fromSequence, sequence: first.fromSequence, index: idx++ });
        seen.add(first.fromSequence);
      }
    }
    for (const s of steps) {
      if (!seen.has(s.toSequence)) {
        nodes.push({ id: s.toSequence, sequence: s.toSequence, index: idx++ });
        seen.add(s.toSequence);
      }
      links.push({ source: s.fromSequence, target: s.toSequence, cost: s.transitionCost });
    }
    return { steps, nodes, links };
  }
}
