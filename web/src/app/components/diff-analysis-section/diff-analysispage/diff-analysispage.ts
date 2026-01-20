import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { DiffAnalysis } from '../diff-analysis/diff-analysis';

@Component({
  selector: 'app-diff-analysispage',
  imports: [
    CommonModule,
    FormsModule,
    DiffAnalysis,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './diff-analysispage.html',
  styleUrl: './diff-analysispage.scss',
  standalone: true
})
export class DiffAnalysispage {}
