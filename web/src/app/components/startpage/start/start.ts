import { Component } from '@angular/core';
import { Count } from '../count/count';
import { CommonModule } from '@angular/common';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';
import { Recount } from '../recount/recount';
import { Preprocess } from '../preprocess/preprocess';

@Component({
  selector: 'app-start',
  imports: [
    CommonModule, 
    Count,
    Recount,
    Preprocess,
    ...MATERIAL_IMPORTS],
  templateUrl: './start.html',
  styleUrl: './start.scss'
})
export class Start {
}
