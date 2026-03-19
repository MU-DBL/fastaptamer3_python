import { Component, PLATFORM_ID, Inject, signal, WritableSignal } from '@angular/core';
import { isPlatformBrowser } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { ClusterDiversity } from '../cluster-diversity/cluster-diversity';
import { ClusterMsa } from '../cluster-msa/cluster-msa';
import { ClusterPhmm } from '../cluster-phmm/cluster-phmm';
import { Cluster } from '../cluster/cluster';
import { Recluster } from '../recluster/recluster';

@Component({
  selector: 'app-clusterpage',
  imports: [
    CommonModule,
    FormsModule,
    Cluster,
    ClusterDiversity,
    ClusterMsa,
    ClusterPhmm,
    Recluster,
    ...MATERIAL_IMPORTS
],
  templateUrl: './clusterpage.html',
  styleUrl: './clusterpage.scss',
})

export class Clusterpage {}
  