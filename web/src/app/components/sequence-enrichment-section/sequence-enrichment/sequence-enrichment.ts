import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

@Component({
  selector: 'app-sequence-enrichment',
  imports: [
    CommonModule,
    FormsModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './sequence-enrichment.html',
  styleUrl: './sequence-enrichment.scss',
  standalone: true
})
export class SequenceEnrichment {
  // TODO: Implement sequence enrichment functionality
  // Based on seqEnrichTab.R from FASTAptameR3
}
