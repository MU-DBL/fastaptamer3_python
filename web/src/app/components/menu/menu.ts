import { Component } from '@angular/core';
import { Router } from '@angular/router';
import { CommonModule } from '@angular/common';

@Component({
  selector: 'app-menu',
  imports: [CommonModule],
  templateUrl: './menu.html',
  styleUrl: './menu.scss'
})

export class Menu {
activeRoute: string = 'Start';

  menuItems = [
    { label: 'Start', route: '/start' },
    { label: 'Translate', route: '/translate' },
    { label: 'Motif', route: '/motif' },
    { label: 'Mutation Network', route: '/mutation-network' },
    { label: 'Distance', route: '/distance' },
    { label: 'Data Merge', route: '/data-merge' },
    { label: 'Sequence Enrichment', route: '/sequence-enrichment' },
    { label: 'Differential Analysis', route: '/differential-analysis' },
    { label: 'Cluster', route: '/cluster' },
    { label: 'About', route: '/about' }
  ];

  constructor(private router: Router) {
    // Set active route based on current URL
    this.activeRoute = this.getCurrentRoute();
  }

  navigateTo(item: any): void {
    this.activeRoute = item.label;
    this.router.navigate([item.route]);
  }

  isActive(label: string): boolean {
    return this.activeRoute === label;
  }

  private getCurrentRoute(): string {
    const currentUrl = this.router.url;
    const menuItem = this.menuItems.find(item => item.route === currentUrl);
    return menuItem ? menuItem.label : 'Start';
  }
}
