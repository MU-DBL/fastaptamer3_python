import { ComponentFixture, TestBed } from '@angular/core/testing';

import { Cluster } from './cluster';

describe('Cluster', () => {
  let component: Cluster;
  let fixture: ComponentFixture<Cluster>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [Cluster]
    })
    .compileComponents();

    fixture = TestBed.createComponent(Cluster);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
