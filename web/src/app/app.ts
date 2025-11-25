import { Component, signal } from '@angular/core';
import { RouterOutlet } from '@angular/router';
import { Menu } from './components/menu/menu';
import { PlotModal } from './components/common/plot-modal/plot-modal';

@Component({
  selector: 'app-root',
  imports: [RouterOutlet, Menu, PlotModal],
  templateUrl: './app.html',
  styleUrl: './app.scss'
})
export class App {
  protected readonly title = signal('Fastaptamer3');
}
