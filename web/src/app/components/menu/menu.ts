import { Component, OnInit } from '@angular/core';
import { Router, NavigationEnd } from '@angular/router';
import { CommonModule } from '@angular/common';
import { filter } from 'rxjs/operators';

@Component({
  selector: 'app-menu',
  imports: [CommonModule],
  templateUrl: './menu.html',
  styleUrl: './menu.scss'
})

export class Menu implements OnInit {
  activeRoute: string = 'Start';
  menuOpen: boolean = false;

  toggleMenu(): void {
    this.menuOpen = !this.menuOpen;
  }

  menuItems = [
    { label: 'Start', route: '/start' },
    { label: 'Translate', route: '/translate' },
    { label: 'Motif', route: '/motif' },
    { label: 'Mutation Network', route: '/mutation-network' },
    { label: 'Distance', route: '/distance' },
    { label: 'Data Merge', route: '/data-merge' },
    { label: 'Enrichment', route: '/enrichment' },
    { label: 'Differential Analysis', route: '/diff-analysis' },
    { label: 'Cluster', route: '/cluster' },
    { label: 'Pipeline', route: '/pipeline' },
    { label: 'About', route: '/about' }
  ];

  constructor(private router: Router) {}

  ngOnInit(): void {
    // Set active route based on current URL
    this.activeRoute = this.getCurrentRoute();

    // Subscribe to router events to update active route on navigation
    this.router.events.pipe(
      filter(event => event instanceof NavigationEnd)
    ).subscribe(() => {
      this.activeRoute = this.getCurrentRoute();
    });
  }

  navigateTo(item: any): void {
    this.activeRoute = item.label;
    this.menuOpen = false;
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
