import { Routes } from '@angular/router';
import { Start }from './components/start-section/start/start'
import { NotFound } from './components/not-found/not-found';
import { Clusterpage } from './components/cluster-section/clusterpage/clusterpage';
import { Motifpage } from './components/motif-section/motifpage/motifpage';
import { Translate } from './components/translate-section/translate/translate';
import { Enrichmentpage } from './components/enrichment-section/enrichmentpage/enrichmentpage';
import { DataMerge } from './components/data-merge-section/data-merge/data-merge';
import { DiffAnalysispage } from './components/diff-analysis-section/diff-analysispage/diff-analysispage';
import { Distance } from './components/distance-section/distance/distance';
import { MutationNetwork } from './components/mutation-network-section/mutation-network/mutation-network';
import { Pipelinepage } from './components/pipeline-section/pipelinepage/pipelinepage';
import { AboutPage } from './components/about/about';

export const routes: Routes = [
  { path: '', component: Start },   // Default route
  { path: 'start', component: Start },    
  { path: 'cluster', component: Clusterpage },
  { path: 'motif', component: Motifpage },
  { path: 'translate', component: Translate },
  { path: 'enrichment', component: Enrichmentpage },
  { path: 'data-merge', component: DataMerge },
  { path: 'diff-analysis', component: DiffAnalysispage },
  { path: 'distance', component: Distance },
  { path: 'mutation-network', component: MutationNetwork },
  { path: 'pipeline', component: Pipelinepage },
  { path: 'about', component: AboutPage },
  { path: '**', component: NotFound }  // 404 - Must be last!
];
