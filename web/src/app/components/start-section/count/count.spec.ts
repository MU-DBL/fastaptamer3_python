import { ComponentFixture, TestBed } from '@angular/core/testing';

import { Count } from './count';

describe('Count', () => {
  let component: Count;
  let fixture: ComponentFixture<Count>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [Count]
    })
    .compileComponents();

    fixture = TestBed.createComponent(Count);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
