import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

@Component({
  selector: 'app-mutation-network',
  imports: [
    CommonModule,
    FormsModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './mutation-network.html',
  styleUrl: './mutation-network.scss',
  standalone: true
})
export class MutationNetwork {
  // TODO: Implement mutation network analysis functionality
  // Based on mutationNetworkTab.R from FASTAptameR3
}
