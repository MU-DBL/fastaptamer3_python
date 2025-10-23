import { Routes } from '@angular/router';
import { Start }from './components/startpage/start/start'
import { NotFound } from './components/not-found/not-found';

export const routes: Routes = [
  { path: '', component: Start },   // Default route
  { path: 'start', component: Start },                            
  { path: '**', component: NotFound }  // 404 - Must be last!
];
