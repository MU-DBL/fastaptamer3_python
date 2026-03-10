import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../shared/material-imports';

@Component({
  selector: 'app-about',
  imports: [CommonModule, ...MATERIAL_IMPORTS],
  templateUrl: './about.html',
  styleUrl: './about.scss'
})
export class AboutPage {}
