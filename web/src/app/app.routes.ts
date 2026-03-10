import { Routes } from '@angular/router';
import { Start }from './components/start-section/start/start'
import { NotFound } from './components/not-found/not-found';
import { Clusterpage } from './components/cluster-section/clusterpage/clusterpage';
import { Motifpage } from './components/motif-section/motifpage/motifpage';
import { Translatepage } from './components/translate-section/translatepage/translatepage';
import { SequenceEnrichmentpage } from './components/sequence-enrichment-section/sequence-enrichmentpage/sequence-enrichmentpage';
import { DataMergepage } from './components/data-merge-section/data-mergepage/data-mergepage';
import { DiffAnalysispage } from './components/diff-analysis-section/diff-analysispage/diff-analysispage';
import { Distancepage } from './components/distance-section/distancepage/distancepage';
import { MutationNetworkpage } from './components/mutation-network-section/mutation-networkpage/mutation-networkpage';
import { Pipelinepage } from './components/pipeline-section/pipelinepage/pipelinepage';
import { AboutPage } from './components/about/about';

export const routes: Routes = [
  { path: '', component: Start },   // Default route
  { path: 'start', component: Start },    
  { path: 'cluster', component: Clusterpage },
  { path: 'motif', component: Motifpage },
  { path: 'translate', component: Translatepage },
  { path: 'sequence-enrichment', component: SequenceEnrichmentpage },
  { path: 'data-merge', component: DataMergepage },
  { path: 'diff-analysis', component: DiffAnalysispage },
  { path: 'distance', component: Distancepage },
  { path: 'mutation-network', component: MutationNetworkpage },
  { path: 'pipeline', component: Pipelinepage },
  { path: 'about', component: AboutPage },
  { path: '**', component: NotFound }  // 404 - Must be last!
];
