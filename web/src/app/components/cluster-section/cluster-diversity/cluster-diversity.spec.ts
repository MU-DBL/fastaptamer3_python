import { ComponentFixture, TestBed } from '@angular/core/testing';

import { ClusterDiversity } from './cluster-diversity';

describe('ClusterDiversity', () => {
  let component: ClusterDiversity;
  let fixture: ComponentFixture<ClusterDiversity>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [ClusterDiversity]
    })
    .compileComponents();

    fixture = TestBed.createComponent(ClusterDiversity);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
