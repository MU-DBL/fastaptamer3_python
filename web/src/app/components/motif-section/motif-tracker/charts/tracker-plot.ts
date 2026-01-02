import { Component, Input, AfterViewInit, OnChanges, SimpleChanges, ElementRef, ViewChild, PLATFORM_ID, inject } from '@angular/core';
import { CommonModule, isPlatformBrowser } from '@angular/common';

interface TrackerDataPoint {
  PopulationNumber: number;
  PopulationName: string;
  [key: string]: any; // Dynamic motif/sequence and percentage/RPU columns
}

@Component({
  selector: 'app-tracker-plot',
  standalone: true,
  imports: [CommonModule],
  template: `
    <div #plotContainer class="chart-container"></div>
  `,
  styles: [`
    .chart-container {
      width: 100%;
      height: 600px;
    }
  `]
})
export class TrackerPlot implements AfterViewInit, OnChanges {
  private readonly platformId = inject(PLATFORM_ID);
  
  @ViewChild('plotContainer', { static: false }) plotContainer!: ElementRef;
  
  @Input() data: TrackerDataPoint[] = [];
  @Input() queryType: 'Motif' | 'Sequence' = 'Motif'; // Determines which column to use for value
  @Input() xAxisLabel = 'Population';
  @Input() yAxisLabel = 'Percentage / RPU';
  @Input() title = 'Motif / Sequence Tracker';
  @Input() colorPalette = 'Dark2'; // Color palette name

  ngAfterViewInit(): void {
    if (isPlatformBrowser(this.platformId)) {
      this.drawChart();
    }
  }

  ngOnChanges(changes: SimpleChanges): void {
    if (isPlatformBrowser(this.platformId) && this.plotContainer) {
      if (changes['data'] || changes['title'] || 
          changes['xAxisLabel'] || changes['yAxisLabel'] || 
          changes['colorPalette'] || changes['queryType']) {
        this.drawChart();
      }
    }
  }

  private async drawChart(): Promise<void> {
    if (!this.plotContainer || !this.data || this.data.length === 0) {
      return;
    }

    try {
      const Plotly = await import('plotly.js-dist-min');

      // Determine which column to use for values (Motif or sequences, depending on tracker type)
      const queryColumn = this.queryType === 'Motif' ? 'Motif' : 'sequences';
      const valueColumn = this.queryType === 'Motif' ? 'Percentage' : 'RPU';
      const aliasColumn = 'Alias';

      // Group data by query (Motif or Sequence)
      const groupedByQuery = this.groupByQuery(this.data, queryColumn);

      // Create color palette
      const colors = this.getColorPalette(this.colorPalette, Object.keys(groupedByQuery).length);

      // Create traces for each query (motif or sequence)
      const traces: any[] = [];
      let colorIndex = 0;

      for (const [query, items] of Object.entries(groupedByQuery)) {
        // Sort by population number to maintain order
        items.sort((a, b) => a.PopulationNumber - b.PopulationNumber);

        // Extract x and y values - use PopulationNumber for proper ordering
        const xValues = items.map(item => item.PopulationNumber);
        const yValues = items.map(item => parseFloat(item[valueColumn]) || 0);
        const populationNames = items.map(item => item.PopulationName);

        // Get alias if available (use first item's alias)
        const alias = items[0]?.[aliasColumn] || query;

        const trace: any = {
          x: xValues,
          y: yValues,
          type: 'scatter',
          mode: 'lines+markers',
          name: alias !== '-' ? alias : query,
          line: {
            color: colors[colorIndex % colors.length],
            width: 3
          },
          marker: {
            size: 8,
            color: colors[colorIndex % colors.length]
          },
          hovertemplate: `<b>${queryColumn}: ${query}</b><br>` +
                        '<b>Population:</b> %{customdata}<br>' +
                        `<b>${valueColumn}:</b> %{y:.2f}<extra></extra>`,
          customdata: populationNames
        };

        traces.push(trace);
        colorIndex++;
      }

      // Create tick labels mapping population numbers to names
      const allPopulations = this.getOrderedPopulations(this.data);
      const tickvals = allPopulations.map((_, idx) => idx + 1); // PopulationNum values
      const ticktext = allPopulations; // PopulationName labels

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
          tickvals: tickvals,
          ticktext: ticktext,
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
        showlegend: true,
        legend: {
          orientation: 'h',
          y: -0.2,
          x: 0.5,
          xanchor: 'center',
          yanchor: 'top'
        },
        margin: { t: 80, b: 120, l: 90, r: 40 },
        plot_bgcolor: 'white',
        paper_bgcolor: 'white'
      };

      const config: any = {
        responsive: true,
        displayModeBar: true,
        displaylogo: false,
        toImageButtonOptions: {
          format: 'svg',
          filename: `${this.queryType.toLowerCase()}_tracker`,
          height: 600,
          width: 900
        },
        modeBarButtonsToAdd: ['pan2d', 'select2d', 'lasso2d'],
        modeBarButtonsToRemove: []
      };

      Plotly.newPlot(this.plotContainer.nativeElement, traces, layout, config);
    } catch (error) {
      console.error('Error creating tracker plot:', error);
    }
  }

  private groupByQuery(data: TrackerDataPoint[], queryColumn: string): Record<string, TrackerDataPoint[]> {
    return data.reduce((acc, item) => {
      const query = item[queryColumn];
      if (!acc[query]) {
        acc[query] = [];
      }
      acc[query].push(item);
      return acc;
    }, {} as Record<string, TrackerDataPoint[]>);
  }

  private getOrderedPopulations(data: TrackerDataPoint[]): string[] {
    const unique = new Map<number, string>();
    data.forEach(item => {
      unique.set(item.PopulationNumber, item.PopulationName);
    });
    return Array.from(unique.entries())
      .sort((a, b) => a[0] - b[0])
      .map(([_, name]) => name);
  }

  private getColorPalette(palette: string, numColors: number): string[] {
    const palettes: Record<string, string[]> = {
      'Dark2': ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666'],
      'Set1': ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#999999', '#a65628', '#f781bf'],
      'Set2': ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'],
      'Paired': ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00'],
      'Pastel1': ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#ffffcc', '#e5d8bd', '#fddaec'],
      'Pastel2': ['#b3e2cd', '#fdcdac', '#f4cae4', '#e2f0d9', '#fff2ae', '#f1e2cc', '#cccccc', '#e0bbe4'],
      'Accent': ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f', '#bf5b17', '#666666']
    };

    const selectedPalette = palettes[palette] || palettes['Dark2'];
    
    // If we need more colors than the palette has, repeat the palette
    if (numColors <= selectedPalette.length) {
      return selectedPalette.slice(0, numColors);
    }

    const repeated: string[] = [];
    for (let i = 0; i < numColors; i++) {
      repeated.push(selectedPalette[i % selectedPalette.length]);
    }
    return repeated;
  }
}
