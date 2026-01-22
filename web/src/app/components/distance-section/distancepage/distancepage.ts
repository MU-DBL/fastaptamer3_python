import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Distance } from '../distance/distance';

@Component({
  selector: 'app-distancepage',
  imports: [
    CommonModule,
    FormsModule,
    Distance,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './distancepage.html',
  styleUrl: './distancepage.scss',
  standalone: true
})
export class Distancepage {}
