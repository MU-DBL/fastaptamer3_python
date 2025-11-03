import { ComponentFixture, TestBed } from '@angular/core/testing';

import { Clusterpage } from './clusterpage';

describe('Clusterpage', () => {
  let component: Clusterpage;
  let fixture: ComponentFixture<Clusterpage>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [Clusterpage]
    })
    .compileComponents();

    fixture = TestBed.createComponent(Clusterpage);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
