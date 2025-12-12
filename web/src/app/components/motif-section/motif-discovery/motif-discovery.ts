import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';

@Component({
  selector: 'app-motif-discovery',
  imports: [ 
    CommonModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './motif-discovery.html',
  styleUrl: './motif-discovery.scss'
})
export class MotifDiscovery {
  // Placeholder component - functionality to be implemented
}
