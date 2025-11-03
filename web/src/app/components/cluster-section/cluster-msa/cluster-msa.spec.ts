import { ComponentFixture, TestBed } from '@angular/core/testing';

import { ClusterMsa } from './cluster-msa';

describe('ClusterMsa', () => {
  let component: ClusterMsa;
  let fixture: ComponentFixture<ClusterMsa>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [ClusterMsa]
    })
    .compileComponents();

    fixture = TestBed.createComponent(ClusterMsa);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
