import { ComponentFixture, TestBed } from '@angular/core/testing';

import { PositionEnrichment } from './position-enrichment';

describe('PositionEnrichment', () => {
  let component: PositionEnrichment;
  let fixture: ComponentFixture<PositionEnrichment>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [PositionEnrichment]
    })
    .compileComponents();

    fixture = TestBed.createComponent(PositionEnrichment);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
