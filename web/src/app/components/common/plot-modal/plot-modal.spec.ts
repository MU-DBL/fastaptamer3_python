import { ComponentFixture, TestBed } from '@angular/core/testing';

import { PlotModal } from './plot-modal';

describe('PlotModal', () => {
  let component: PlotModal;
  let fixture: ComponentFixture<PlotModal>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [PlotModal]
    })
    .compileComponents();

    fixture = TestBed.createComponent(PlotModal);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
