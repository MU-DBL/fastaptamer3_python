declare module 'plotly.js-dist-min' {
  export function newPlot(
    root: HTMLElement,
    data: any[],
    layout?: any,
    config?: any
  ): Promise<any>;

  export function react(
    root: HTMLElement,
    data: any[],
    layout?: any,
    config?: any
  ): Promise<any>;

  export function relayout(root: HTMLElement, layout: any): Promise<any>;
  export function restyle(root: HTMLElement, update: any, traces?: number | number[]): Promise<any>;
  export function update(root: HTMLElement, dataUpdate: any, layoutUpdate: any, traces?: number | number[]): Promise<any>;
  export function purge(root: HTMLElement): void;
}
