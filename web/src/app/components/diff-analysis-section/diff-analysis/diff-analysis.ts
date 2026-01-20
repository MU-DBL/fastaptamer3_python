import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

@Component({
  selector: 'app-diff-analysis',
  imports: [
    CommonModule,
    FormsModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './diff-analysis.html',
  styleUrl: './diff-analysis.scss',
  standalone: true
})
export class DiffAnalysis {
  // TODO: Implement differential analysis functionality
  // Based on diffAnalysisTab.R from FASTAptameR3
}
