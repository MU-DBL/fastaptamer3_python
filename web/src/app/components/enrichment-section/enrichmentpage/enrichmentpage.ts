import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { SequenceEnrichment } from '../sequence-enrichment/sequence-enrichment';
import { PositionEnrichment } from '../position-enrichment/position-enrichment';

@Component({
  selector: 'app-enrichmentpage',
  imports: [
    CommonModule,
    FormsModule,
    SequenceEnrichment,
    PositionEnrichment,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './enrichmentpage.html',
  styleUrl: './enrichmentpage.scss',
})
export class Enrichmentpage {}
