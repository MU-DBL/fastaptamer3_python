import { ComponentFixture, TestBed } from '@angular/core/testing';

import { ClusterPhmm } from './cluster-phmm';

describe('ClusterPhmm', () => {
  let component: ClusterPhmm;
  let fixture: ComponentFixture<ClusterPhmm>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [ClusterPhmm]
    })
    .compileComponents();

    fixture = TestBed.createComponent(ClusterPhmm);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
