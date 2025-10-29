import { Component, Input, AfterViewInit, OnChanges, SimpleChanges, ElementRef, ViewChild, PLATFORM_ID, inject } from '@angular/core';
import { CommonModule, isPlatformBrowser } from '@angular/common';

interface RgbColor {
  r: number;
  g: number;
  b: number;
}

interface BinnedData {
  bin: string;
  fraction: number;
  uniqueCount: number;
}

@Component({
  selector: 'app-abundance-plot',
  standalone: true,
  imports: [CommonModule],
  template: `
    <div #plotContainer class="chart-container"></div>
  `,
  styles: [`
    .chart-container {
      width: 100%;
      height: 550px;
    }
  `]
})
export class AbundancePlot implements AfterViewInit, OnChanges {
  private readonly platformId = inject(PLATFORM_ID);
  
  @ViewChild('plotContainer', { static: false }) plotContainer!: ElementRef;
  
  @Input() data: BinnedData[] = [];
  @Input() xAxisLabel = 'No. Reads';
  @Input() yAxisLabel = 'Fraction of Population';
  @Input() title = 'Binned sequence abundance';
  @Input() colorLight = '#ADD8E6'; // Light blue for low values
  @Input() colorDark = '#FF6B6B';  // Red/coral for high values

  ngAfterViewInit(): void {
    if (isPlatformBrowser(this.platformId)) {
      this.drawChart();
    }
  }

  ngOnChanges(changes: SimpleChanges): void {
    if (isPlatformBrowser(this.platformId) && this.plotContainer) {
      // Only redraw if data or visual properties changed
      if (changes['data'] || changes['colorLight'] || changes['colorDark'] || 
          changes['title'] || changes['xAxisLabel'] || changes['yAxisLabel']) {
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

      const xValues = this.data.map(d => d.bin);
      const yValues = this.data.map(d => d.fraction);
      const uniqueCounts = this.data.map(d => d.uniqueCount);
      const maxUniqueCount = Math.max(...uniqueCounts);

      const lightRgb = this.hexToRgb(this.colorLight);
      const darkRgb = this.hexToRgb(this.colorDark);

      // Create color array based on unique count intensity
      const colors = this.generateGradientColors(uniqueCounts, maxUniqueCount, lightRgb, darkRgb);

      const trace: any = {
        x: xValues,
        y: yValues,
        type: 'bar',
        marker: {
          color: colors,
          line: {
            color: '#000',
            width: 1
          }
        },
        text: uniqueCounts.map(count => `Unique: ${count}`),
        hovertemplate: '<b>Bin:</b> %{x}<br><b>Fraction:</b> %{y:.4f}<br><b>%{text}</b><extra></extra>'
      };

      const { legendAnnotations, legendShapes } = this.createGradientLegend(maxUniqueCount, lightRgb, darkRgb);

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
          tickangle: -45,
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
        margin: { t: 60, b: 140, l: 90, r: 180 },
        plot_bgcolor: 'white',
        paper_bgcolor: 'white',
        annotations: legendAnnotations,
        shapes: legendShapes
      };

      const config: any = {
        responsive: true,
        displayModeBar: true,
        displaylogo: false,
        toImageButtonOptions: {
          format: 'svg',
          filename: 'binned_abundance',
          height: 500,
          width: 900
        },
        modeBarButtonsToAdd: ['pan2d'],
        modeBarButtonsToRemove: []
      };

      Plotly.newPlot(this.plotContainer.nativeElement, [trace], layout, config);
    } catch (error) {
      console.error('Error creating abundance plot:', error);
    }
  }

  /**
   * Parse hex color string to RGB object
   */
  private hexToRgb(hex: string): RgbColor {
    const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
      r: parseInt(result[1], 16),
      g: parseInt(result[2], 16),
      b: parseInt(result[3], 16)
    } : { r: 173, g: 216, b: 230 }; // Default to light blue
  }

  /**
   * Generate gradient colors interpolating between light and dark colors
   */
  private generateGradientColors(
    uniqueCounts: number[], 
    maxUniqueCount: number, 
    lightRgb: RgbColor, 
    darkRgb: RgbColor
  ): string[] {
    return uniqueCounts.map(count => {
      const intensity = maxUniqueCount > 0 ? count / maxUniqueCount : 0;
      return this.interpolateColor(lightRgb, darkRgb, intensity);
    });
  }

  /**
   * Interpolate between two RGB colors based on intensity (0-1)
   */
  private interpolateColor(lightRgb: RgbColor, darkRgb: RgbColor, intensity: number): string {
    const r = Math.round(lightRgb.r + (darkRgb.r - lightRgb.r) * intensity);
    const g = Math.round(lightRgb.g + (darkRgb.g - lightRgb.g) * intensity);
    const b = Math.round(lightRgb.b + (darkRgb.b - lightRgb.b) * intensity);
    return `rgb(${r}, ${g}, ${b})`;
  }

  /**
   * Create gradient legend with continuous color bar
   */
  private createGradientLegend(
    maxUniqueCount: number, 
    lightRgb: RgbColor, 
    darkRgb: RgbColor
  ): { legendAnnotations: any[], legendShapes: any[] } {
    const legendSteps = 20;
    const legendAnnotations: any[] = [];
    const legendShapes: any[] = [];
    
    // Legend configuration
    const gradientHeight = 0.6;
    const gradientWidth = 0.035;
    const gradientStartY = 0.15;
    
    // Add title for legend
    legendAnnotations.push({
      text: `<b>No. Unique</b><br><b>Sequences</b>`,
      xref: 'paper',
      yref: 'paper',
      x: 1.095,
      y: 0.88,
      xanchor: 'center',
      yanchor: 'bottom',
      showarrow: false,
      font: { size: 10, family: 'Arial' },
      align: 'center'
    });
    
    // Create continuous gradient bar using rectangles
    for (let i = 0; i < legendSteps; i++) {
      const intensity = 1 - (i / (legendSteps - 1)); // Top is darkest
      const color = this.interpolateColor(lightRgb, darkRgb, intensity);
      
      const yStart = gradientStartY + (i / legendSteps) * gradientHeight;
      const yEnd = gradientStartY + ((i + 1) / legendSteps) * gradientHeight;
      
      legendShapes.push({
        type: 'rect',
        xref: 'paper',
        yref: 'paper',
        x0: 1.04,
        y0: yStart,
        x1: 1.04 + gradientWidth,
        y1: yEnd,
        fillcolor: color,
        line: { width: 0 }
      });
    }
    
    // Add value labels at key points
    const labelPositions = [0, 0.25, 0.5, 0.75, 1.0];
    labelPositions.forEach(pos => {
      const intensity = 1 - pos;
      const value = Math.round(maxUniqueCount * intensity);
      const yPos = gradientStartY + pos * gradientHeight;
      
      legendAnnotations.push({
        xref: 'paper',
        yref: 'paper',
        x: 1.12,
        y: yPos,
        xanchor: 'left',
        yanchor: 'middle',
        text: value.toString(),
        showarrow: false,
        font: { size: 9, family: 'Arial' },
        align: 'left'
      });
    });

    return { legendAnnotations, legendShapes };
  }
}
