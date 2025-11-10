/**
 * Color Service for Cluster Visualization
 * 
 * Handles color assignment and palette management for cluster visualization.
 * Follows the RColorBrewer palette standards used in the R implementation.
 */

export interface ClusterColor {
  cluster: number;
  color: string;
}

export interface ColorPalette {
  name: string;
  colors: string[];
  description: string;
}

export class ClusterColorService {
  /**
   * RColorBrewer-inspired color palettes matching R implementation
   */
  private static readonly PALETTES: { [key: string]: ColorPalette } = {
    'Dark2': {
      name: 'Dark2',
      colors: ['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666'],
      description: 'Dark colors for high contrast visualization'
    },
    'Set1': {
      name: 'Set1',
      colors: ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF'],
      description: 'Bright, distinct colors for cluster differentiation'
    },
    'Set2': {
      name: 'Set2',
      colors: ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3'],
      description: 'Softer colors for gentle visualization'
    },
    'Paired': {
      name: 'Paired',
      colors: ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00'],
      description: 'Paired colors for related cluster grouping'
    },
    'Accent': {
      name: 'Accent',
      colors: ['#7FC97F', '#BEAED4', '#FDC086', '#FFFF99', '#386CB0', '#F0027F', '#BF5B17', '#666666'],
      description: 'Accent colors for highlighting specific clusters'
    }
  };

  /**
   * Default color for non-clustered sequences (cluster 0)
   */
  private static readonly NON_CLUSTERED_COLOR = '#999999';

  /**
   * Get available color palettes
   */
  static getAvailablePalettes(): ColorPalette[] {
    return Object.values(this.PALETTES);
  }

  /**
   * Get a specific color palette by name
   */
  static getPalette(paletteName: string): ColorPalette | null {
    return this.PALETTES[paletteName] || null;
  }

  /**
   * Assign colors to clusters using specified palette
   */
  static assignClusterColors(
    clusterIds: number[], 
    paletteName: string = 'Dark2'
  ): ClusterColor[] {
    const palette = this.getPalette(paletteName);
    if (!palette) {
      throw new Error(`Unknown color palette: ${paletteName}`);
    }

    const colors = palette.colors;
    const clusterColors: ClusterColor[] = [];

    clusterIds.forEach((clusterId, index) => {
      let color: string;
      
      if (clusterId === 0) {
        // Non-clustered sequences get a special gray color
        color = this.NON_CLUSTERED_COLOR;
      } else {
        // Assign colors cyclically from the palette
        color = colors[index % colors.length];
      }

      clusterColors.push({
        cluster: clusterId,
        color: color
      });
    });

    return clusterColors;
  }

  /**
   * Create a color mapping object for easy lookup
   */
  static createColorMap(
    clusterIds: number[], 
    paletteName: string = 'Dark2'
  ): { [cluster: number]: string } {
    const clusterColors = this.assignClusterColors(clusterIds, paletteName);
    const colorMap: { [cluster: number]: string } = {};
    
    clusterColors.forEach(({ cluster, color }) => {
      colorMap[cluster] = color;
    });

    return colorMap;
  }

  /**
   * Get color for a specific cluster using the palette
   */
  static getClusterColor(
    clusterId: number, 
    allClusterIds: number[], 
    paletteName: string = 'Dark2'
  ): string {
    if (clusterId === 0) {
      return this.NON_CLUSTERED_COLOR;
    }

    const palette = this.getPalette(paletteName);
    if (!palette) {
      return this.NON_CLUSTERED_COLOR;
    }

    // Find the index of this cluster in the sorted list
    const sortedClusters = allClusterIds.filter(id => id !== 0).sort((a, b) => a - b);
    const index = sortedClusters.indexOf(clusterId);
    
    if (index === -1) {
      return this.NON_CLUSTERED_COLOR;
    }

    return palette.colors[index % palette.colors.length];
  }

  /**
   * Validate palette name
   */
  static isValidPalette(paletteName: string): boolean {
    return paletteName in this.PALETTES;
  }

  /**
   * Get default palette name
   */
  static getDefaultPalette(): string {
    return 'Dark2';
  }
}