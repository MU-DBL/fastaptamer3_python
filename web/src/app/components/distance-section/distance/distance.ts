import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

@Component({
  selector: 'app-distance',
  imports: [
    CommonModule,
    FormsModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './distance.html',
  styleUrl: './distance.scss',
  standalone: true
})
export class Distance {
  // TODO: Implement distance calculation functionality
  // Based on distanceTab.R from FASTAptameR3
}
