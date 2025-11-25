import { Component, OnInit, OnDestroy, ViewChild, ElementRef, ChangeDetectorRef, AfterViewInit, PLATFORM_ID, Inject, NgZone } from '@angular/core';
import { CommonModule, isPlatformBrowser } from '@angular/common';
import { combineLatest, distinctUntilChanged, filter, Subject, takeUntil } from 'rxjs';
import { PlotModalService, PlotConfig } from '../../../shared/plot-modal.service';
import { FormsModule } from '@angular/forms';
import { MATERIAL_IMPORTS } from '../../../shared/material-imports';

@Component({
  selector: 'app-plot-modal',
  standalone: true,
  imports: [
    CommonModule,
    FormsModule,
    ...MATERIAL_IMPORTS
  ],
  templateUrl: './plot-modal.html',
  styleUrls: ['./plot-modal.scss']
})
export class PlotModal implements AfterViewInit, OnDestroy {
  @ViewChild('plotContainer') plotContainer!: ElementRef;

  isOpen = false;
  plotConfig: PlotConfig | null = null;
  private destroy$ = new Subject<void>();
  private PlotlyJS: any;

  constructor(
    private cdr: ChangeDetectorRef,
    private plotModalService: PlotModalService,
    private ngZone: NgZone,
    @Inject(PLATFORM_ID) private platformId: Object
  ) {}

  ngAfterViewInit(): void {
    // ✅ SINGLE subscription that waits for BOTH conditions
    combineLatest([
      this.plotModalService.isOpen$,
      this.plotModalService.plotConfig$
    ])
      .pipe(
        takeUntil(this.destroy$),
        distinctUntilChanged((prev, curr) => 
          prev[0] === curr[0] && prev[1] === curr[1]
        ),
        filter(([isOpen, config]) => isOpen && config !== null)
      )
      .subscribe(([isOpen, config]) => {
        this.isOpen = isOpen;
        this.plotConfig = config;
        this.cdr.detectChanges();
        
        // ✅ Create plot outside Angular zone
        this.ngZone.runOutsideAngular(() => {
          // Microtask ensures DOM is fully updated
          Promise.resolve().then(() => this.createPlot());
        });
      });
  }

  ngOnDestroy(): void {
    if (isPlatformBrowser(this.platformId)) {
        window.removeEventListener('resize', this.handleResize);
    }   

    if (this.plotContainer && this.PlotlyJS) {
      try {
        this.PlotlyJS.purge(this.plotContainer.nativeElement);
      } catch (e) {
        console.warn('Plotly cleanup error:', e);
      }
    }

    this.destroy$.next();
    this.destroy$.complete();
  }

  close(): void {
    if (isPlatformBrowser(this.platformId)) {
        window.removeEventListener('resize', this.handleResize);
    }
    if (this.plotContainer && this.PlotlyJS) {
      try {
          this.PlotlyJS.purge(this.plotContainer.nativeElement);
      } catch (e) {
          console.warn('Purge error', e);
      }
    }
    
    this.plotModalService.close();
    this.isOpen = false;
  }

  private async createPlot(): Promise<void> {
    if (!isPlatformBrowser(this.platformId)) return;
    if (!this.plotConfig || !this.plotContainer?.nativeElement) {
      console.error('Cannot create plot: missing config or container');
      return;
    }

    if (!this.PlotlyJS) {
      try {
        this.PlotlyJS = await import('plotly.js-dist-min');
      } catch (error) {
        console.error('Failed to load Plotly:', error);
        return;
      }
    }

    const container = this.plotContainer.nativeElement;
    
    try {
      this.PlotlyJS.purge(container);
    } catch (e) {
      // First render, no plot to purge
    }

    if (container.offsetWidth === 0) {
      console.warn('Container not yet visible, retrying...');
      setTimeout(() => this.createPlot(), 50);
      return;
    }

    const containerWidth = container.offsetWidth;
    const containerHeight = container.offsetHeight;

    console.log('Creating plot:', containerWidth, 'x', containerHeight);

    const layout = {
      ...this.plotConfig.layout,
      autosize: true,
      width: containerWidth,
      height: containerHeight,
      margin: { l: 60, r: 40, t: 60, b: 60 },
    };

    const defaultConfig = { 
      responsive: true,
      displayModeBar: true,
      displaylogo: false,
      modeBarButtonsToRemove: ['lasso2d', 'select2d'],
      toImageButtonOptions: { 
        format: 'svg' as const, 
        filename: this.getPlotTitle(),
        width: containerWidth,
        height: containerHeight
      }
    };

    const config = { ...defaultConfig, ...this.plotConfig.config };

    try {
      await this.PlotlyJS.newPlot(container, this.plotConfig.data, layout, config);
      console.log('Plot created successfully');

      window.removeEventListener('resize', this.handleResize);
      window.addEventListener('resize', this.handleResize);
      
    } catch (err: any) {
      console.error('Plotly error:', err);
    }
  }

  private handleResize = (): void => {
    if (!this.plotContainer?.nativeElement || !this.PlotlyJS) return;
    
    try {
      this.PlotlyJS.Plots.resize(this.plotContainer.nativeElement);
    } catch (e) {
      console.warn('Resize error:', e);
    }
  };

  getPlotTitle(): string {
    const title = this.plotConfig?.layout?.title;
    if (typeof title === 'string') return title;
    return (title as any)?.text || 'plot';
  }
}