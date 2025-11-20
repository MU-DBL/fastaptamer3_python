import { Component, Input, AfterViewInit, OnChanges, SimpleChanges, ElementRef, ViewChild, PLATFORM_ID, inject } from '@angular/core';
import { CommonModule, isPlatformBrowser } from '@angular/common';

interface SequenceLengthData {
  unique: { length: number; count: number }[];
  total: { length: number; count: number }[];
}

@Component({
  selector: 'app-seq-length-histogram',
  standalone: true,
  imports: [CommonModule],
  template: `
    <div #plotContainer class="chart-container"></div>
  `,
  styles: [`
    .chart-container {
      width: 100%;
      height: 700px;
    }
  `]
})
export class SeqLengthHistogram implements AfterViewInit, OnChanges {
  private readonly platformId = inject(PLATFORM_ID);
  
  @ViewChild('plotContainer', { static: false }) plotContainer!: ElementRef;
  
  @Input() data: SequenceLengthData = { unique: [], total: [] };
  @Input() xAxisLabel = 'Sequence length (nucleotides)';
  @Input() yAxis1Label = 'Number of unique sequences';
  @Input() yAxis2Label = 'Total reads';
  @Input() title = 'Distribution of Sequence Lengths';
  @Input() uniqueColor = '#87CEEB';
  @Input() totalColor = '#FFA500';

  ngAfterViewInit(): void {
    if (isPlatformBrowser(this.platformId)) {
      this.drawChart();
    }
  }

  ngOnChanges(changes: SimpleChanges): void {
    if (isPlatformBrowser(this.platformId) && this.plotContainer) {
      // Only redraw if data or visual properties changed
      if (changes['data'] || changes['uniqueColor'] || changes['totalColor'] || 
          changes['title'] || changes['xAxisLabel'] || changes['yAxis1Label'] || changes['yAxis2Label']) {
        this.drawChart();
      }
    }
  }

  private async drawChart(): Promise<void> {
    if (!this.plotContainer || !this.data.unique || this.data.unique.length === 0) {
      return;
    }

    try {
      // Dynamically import Plotly only in the browser
      const Plotly = await import('plotly.js-dist-min');

      const lengths = this.data.unique.map(d => d.length);
      const uniqueValues = this.data.unique.map(d => d.count);
      const totalValues = this.data.total.map(d => d.count);

      // Trace for unique sequences
      const traceUnique: any = {
        x: lengths,
        y: uniqueValues,
        type: 'bar',
        name: this.yAxis1Label,
        marker: {
          color: this.uniqueColor,
          line: {
            color: '#000',
            width: 0.5
          }
        },
        hovertemplate: '<b>Length:</b> %{x}<br><b>Unique:</b> %{y}<extra></extra>'
      };

      // Trace for total reads
      const traceTotal: any = {
        x: lengths,
        y: totalValues,
        type: 'bar',
        name: this.yAxis2Label,
        marker: {
          color: this.totalColor,
          line: {
            color: '#000',
            width: 0.5
          }
        },
        hovertemplate: '<b>Length:</b> %{x}<br><b>Total:</b> %{y}<extra></extra>',
        xaxis: 'x2',
        yaxis: 'y2'
      };

      const layout: any = {
        title: {
          text: `<b>${this.title}</b>`,
          font: { size: 18 }
        },
        grid: {
          rows: 2,
          columns: 1,
          pattern: 'independent',
          roworder: 'top to bottom'
        },
        xaxis: {
          title: '',
          showgrid: true,
          gridcolor: '#e0e0e0',
          showline: true,
          linewidth: 2,
          linecolor: 'black',
          mirror: true,
          domain: [0, 1],
          anchor: 'y'
        },
        yaxis: {
          title: {
            text: `<b>${this.yAxis1Label}</b>`,
            font: { size: 14 }
          },
          showgrid: true,
          gridcolor: '#e0e0e0',
          showline: true,
          linewidth: 2,
          linecolor: 'black',
          mirror: true,
          domain: [0.53, 1],
          anchor: 'x'
        },
        xaxis2: {
          title: {
            text: `<b>${this.xAxisLabel}</b>`,
            font: { size: 14 }
          },
          showgrid: true,
          gridcolor: '#e0e0e0',
          showline: true,
          linewidth: 2,
          linecolor: 'black',
          mirror: true,
          domain: [0, 1],
          anchor: 'y2'
        },
        yaxis2: {
          title: {
            text: `<b>${this.yAxis2Label}</b>`,
            font: { size: 14 }
          },
          showgrid: true,
          gridcolor: '#e0e0e0',
          showline: true,
          linewidth: 2,
          linecolor: 'black',
          mirror: true,
          domain: [0, 0.47],
          anchor: 'x2'
        },
        showlegend: false,
        hovermode: 'closest',
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
          filename: 'sequence_length_histogram',
          height: 700,
          width: 900
        },
        modeBarButtonsToAdd: ['pan2d'],
        modeBarButtonsToRemove: []
      };

      Plotly.newPlot(this.plotContainer.nativeElement, [traceUnique, traceTotal], layout, config);
    } catch (error) {
      console.error('Error creating sequence length histogram:', error);
    }
  }
}
