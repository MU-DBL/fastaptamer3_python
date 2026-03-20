import { Component, ElementRef, HostListener, Input, ViewChild } from '@angular/core';
import { CommonModule } from '@angular/common';

@Component({
  selector: 'app-split-panel',
  imports: [CommonModule],
  templateUrl: './split-panel.html',
  styleUrl: './split-panel.scss'
})
export class SplitPanel {
  @Input() initialLeftPercent: number = 30;
  @ViewChild('container') containerRef!: ElementRef<HTMLElement>;

  leftPercent: number = this.initialLeftPercent;
  isDragging = false;

  ngOnInit(): void {
    this.leftPercent = this.initialLeftPercent;
  }

  onDividerMouseDown(event: MouseEvent): void {
    event.preventDefault();
    this.isDragging = true;
  }

  @HostListener('document:mousemove', ['$event'])
  onMouseMove(event: MouseEvent): void {
    if (!this.isDragging) return;
    const container = this.containerRef.nativeElement;
    const rect = container.getBoundingClientRect();
    const offsetX = event.clientX - rect.left;
    const newPercent = (offsetX / rect.width) * 100;
    this.leftPercent = Math.min(Math.max(newPercent, 15), 75);
  }

  @HostListener('document:mouseup')
  onMouseUp(): void {
    this.isDragging = false;
  }
}
