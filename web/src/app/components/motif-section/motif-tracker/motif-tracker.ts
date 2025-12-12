import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';

@Component({
  selector: 'app-motif-tracker',
  imports: [ 
    CommonModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './motif-tracker.html',
  styleUrl: './motif-tracker.scss'
})
export class MotifTracker {
  // Placeholder component - functionality to be implemented
}
