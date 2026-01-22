import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { SequenceEnrichment } from '../sequence-enrichment/sequence-enrichment';

@Component({
  selector: 'app-sequence-enrichmentpage',
  imports: [
    CommonModule,
    FormsModule,
    SequenceEnrichment,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './sequence-enrichmentpage.html',
  styleUrl: './sequence-enrichmentpage.scss',
  standalone: true
})
export class SequenceEnrichmentpage {}
