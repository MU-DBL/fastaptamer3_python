// services/plot-modal.service.ts
import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable } from 'rxjs';

export interface PlotConfig {
  data: any[];
  layout: any;  // Plotly layout object
  config?: any; // Optional Plotly config
}

@Injectable({
  providedIn: 'root'
})
export class PlotModalService {
  private plotConfigSubject = new BehaviorSubject<PlotConfig | null>(null);
  private isOpenSubject = new BehaviorSubject<boolean>(false);

  public plotConfig$: Observable<PlotConfig | null> = this.plotConfigSubject.asObservable();
  public isOpen$: Observable<boolean> = this.isOpenSubject.asObservable();

  openPlot(config: PlotConfig): void {
    console.log('Opening plot:', config.layout?.title);
    this.plotConfigSubject.next(config);
    this.isOpenSubject.next(true);
  }

  close(): void {
    this.isOpenSubject.next(false);
    // Remove delayed cleanup to fix double-click issue
    this.plotConfigSubject.next(null);
  }

  isOpen(): boolean {
    return this.isOpenSubject.value;
  }
}