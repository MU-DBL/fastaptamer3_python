import { Component, Input, AfterViewInit, OnChanges, SimpleChanges, ElementRef, ViewChild, PLATFORM_ID, inject } from '@angular/core';
import { CommonModule, isPlatformBrowser } from '@angular/common';

interface RankData {
  rank: number;
  reads: number;
}

@Component({
  selector: 'app-reads-per-rank-chart',
  standalone: true,
  imports: [CommonModule],
  template: `
    <div #plotContainer class="chart-container"></div>
  `,
  styles: [`
    .chart-container {
      width: 100%;
      height: 500px;
    }
  `]
})
export class ReadsPerRankChart implements AfterViewInit, OnChanges {
  private readonly platformId = inject(PLATFORM_ID);
  
  @ViewChild('plotContainer', { static: false }) plotContainer!: ElementRef;
  
  @Input() data: RankData[] = [];
  @Input() xAxisLabel = 'Ranks of unique sequences';
  @Input() yAxisLabel = 'Total reads per unique sequence';
  @Input() title = 'Read count for each rank';
  @Input() lineColor = '#87CEEB';

  ngAfterViewInit(): void {
    if (isPlatformBrowser(this.platformId)) {
      this.drawChart();
    }
  }

  ngOnChanges(changes: SimpleChanges): void {
    if (isPlatformBrowser(this.platformId) && this.plotContainer) {
      // Only redraw if data or visual properties changed
      if (changes['data'] || changes['lineColor'] || changes['title'] || 
          changes['xAxisLabel'] || changes['yAxisLabel']) {
        this.drawChart();
      }
    }
  }

  private async drawChart(): Promise<void> {
    if (!this.plotContainer || !this.data || this.data.length === 0) {
      return;
    }

    try {
      // Dynamically import Plotly only in the browser
      const Plotly = await import('plotly.js-dist-min');

      const xValues = this.data.map(d => d.rank);
      const yValues = this.data.map(d => d.reads);

      const trace: any = {
        x: xValues,
        y: yValues,
        type: 'scatter',
        mode: 'lines',
        line: {
          color: this.lineColor,
          width: 3
        },
        hovertemplate: '<b>Rank:</b> %{x}<br><b>Reads:</b> %{y}<extra></extra>'
      };

      const layout: any = {
        title: {
          text: `<b>${this.title}</b>`,
          font: { size: 18 }
        },
        xaxis: {
          title: {
            text: `<b>${this.xAxisLabel}</b>`,
            font: { size: 14 }
          },
          showgrid: true,
          gridcolor: '#e0e0e0',
          showline: true,
          linewidth: 2,
          linecolor: 'black',
          mirror: true
        },
        yaxis: {
          title: {
            text: `<b>${this.yAxisLabel}</b>`,
            font: { size: 14 }
          },
          showgrid: true,
          gridcolor: '#e0e0e0',
          showline: true,
          linewidth: 2,
          linecolor: 'black',
          mirror: true
        },
        hovermode: 'closest',
        showlegend: false,
        margin: { t: 60, b: 70, l: 90, r: 40 },
        plot_bgcolor: 'white',
        paper_bgcolor: 'white'
      };

      const config: any = {
        responsive: true,
        displayModeBar: true,
        displaylogo: false,
        toImageButtonOptions: {
          format: 'svg',
          filename: 'reads_per_rank',
          height: 500,
          width: 900
        },
        modeBarButtonsToAdd: ['pan2d', 'select2d', 'lasso2d'],
        modeBarButtonsToRemove: []
      };

      Plotly.newPlot(this.plotContainer.nativeElement, [trace], layout, config);
    } catch (error) {
      console.error('Error creating reads per rank chart:', error);
    }
  }
}
