import { Component } from '@angular/core';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Translate } from '../translate/translate';

@Component({
  selector: 'app-translatepage',
  imports: [
    CommonModule,
    FormsModule,
    Translate,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './translatepage.html',
  styleUrl: './translatepage.scss',
})
export class Translatepage {}
