import { Component, signal, OnInit, PLATFORM_ID, inject } from '@angular/core';
import { isPlatformBrowser } from '@angular/common';
import { RouterOutlet } from '@angular/router';
import { Menu } from './components/menu/menu';
import { PlotModal } from './components/common/plot-modal/plot-modal';
import { environment } from '../environments/environment';

@Component({
  selector: 'app-root',
  imports: [RouterOutlet, Menu, PlotModal],
  templateUrl: './app.html',
  styleUrl: './app.scss'
})
export class App implements OnInit {
  protected readonly title = signal('Fastaptamer3');
  private readonly platformId = inject(PLATFORM_ID);

  ngOnInit(): void {
    if (isPlatformBrowser(this.platformId)) {
      // Clear any leftover files from the previous session on startup
      fetch(`${environment.apiUrl}/clear-files`, { method: 'POST' });
    }
  }
}
