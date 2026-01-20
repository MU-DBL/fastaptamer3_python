import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MutationNetwork } from '../mutation-network/mutation-network';

@Component({
  selector: 'app-mutation-networkpage',
  imports: [
    CommonModule,
    FormsModule,
    MutationNetwork,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './mutation-networkpage.html',
  styleUrl: './mutation-networkpage.scss',
  standalone: true
})
export class MutationNetworkpage {}
