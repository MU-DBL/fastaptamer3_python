import { ComponentFixture, TestBed } from '@angular/core/testing';

import { Recluster } from './recluster';

describe('Recluster', () => {
  let component: Recluster;
  let fixture: ComponentFixture<Recluster>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [Recluster]
    })
    .compileComponents();

    fixture = TestBed.createComponent(Recluster);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
