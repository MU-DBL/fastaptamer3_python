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

  private Plotly: any;

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
      const { default: Plotly } = await import('plotly.js-dist-min');
      this.Plotly = Plotly;

      // Sort by rank and break the line wherever ranks aren't consecutive,
      // so filtered-out/omitted ranks show as a gap instead of a misleading
      // straight line connecting the surrounding points.
      const sorted = [...this.data].sort((a, b) => a.rank - b.rank);
      const xValues: (number | null)[] = [];
      const yValues: (number | null)[] = [];

      // Rank is always 1-indexed by definition, so a missing rank 1 is a
      // leading gap with no preceding point to break against. Anchor the
      // trace (and therefore the x-axis range) at rank 1 with a null value
      // so the missing leading ranks are visible instead of the axis just
      // autoranging to start at the first surviving rank.
      if (sorted.length > 0 && sorted[0].rank > 1) {
        xValues.push(1);
        yValues.push(null);
      }

      for (let i = 0; i < sorted.length; i++) {
        xValues.push(sorted[i].rank);
        yValues.push(sorted[i].reads);
        if (i < sorted.length - 1 && sorted[i + 1].rank - sorted[i].rank > 1) {
          xValues.push(sorted[i].rank + 1);
          yValues.push(null);
        }
      }

      // With gaps now breaking the line, a real point whose neighbors on
      // both sides are broken (e.g. rank 1 survives but rank 2 doesn't) has
      // no line segment touching it and would be invisible in mode:'lines'.
      // Give only those isolated points a visible marker so they still show.
      const markerSizes = yValues.map((y, i) => {
        if (y === null) return 0;
        const leftBroken = i === 0 || yValues[i - 1] === null;
        const rightBroken = i === yValues.length - 1 || yValues[i + 1] === null;
        return leftBroken && rightBroken ? 6 : 0;
      });

      const trace: any = {
        x: xValues,
        y: yValues,
        type: 'scatter',
        mode: 'lines+markers',
        connectgaps: false,
        marker: {
          size: markerSizes,
          color: this.lineColor
        },
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
