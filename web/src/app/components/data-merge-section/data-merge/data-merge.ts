import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

@Component({
  selector: 'app-data-merge',
  imports: [
    CommonModule,
    FormsModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './data-merge.html',
  styleUrl: './data-merge.scss',
  standalone: true
})
export class DataMerge {
  // TODO: Implement data merge functionality
  // Based on dataMergeTab.R from FASTAptameR3
}
