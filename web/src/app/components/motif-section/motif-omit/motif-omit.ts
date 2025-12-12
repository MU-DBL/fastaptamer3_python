import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';

@Component({
  selector: 'app-motif-omit',
  imports: [ 
    CommonModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './motif-omit.html',
  styleUrl: './motif-omit.scss'
})
export class MotifOmit {
  // Placeholder component - functionality to be implemented
}
