import { Routes } from '@angular/router';
import { Start }from './components/start-section/start/start'
import { NotFound } from './components/not-found/not-found';
import { Clusterpage } from './components/cluster-section/clusterpage/clusterpage';

export const routes: Routes = [
  { path: '', component: Start },   // Default route
  { path: 'start', component: Start },    
  { path: 'cluster', component: Clusterpage },                                  
  { path: '**', component: NotFound }  // 404 - Must be last!
];
