import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { DataMerge } from '../data-merge/data-merge';

@Component({
  selector: 'app-data-mergepage',
  imports: [
    CommonModule,
    FormsModule,
    DataMerge,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './data-mergepage.html',
  styleUrl: './data-mergepage.scss',
  standalone: true
})
export class DataMergepage {}
