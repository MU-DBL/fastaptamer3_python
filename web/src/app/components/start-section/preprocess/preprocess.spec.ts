import { ComponentFixture, TestBed } from '@angular/core/testing';

import { Preprocess } from './preprocess';

describe('Preprocess', () => {
  let component: Preprocess;
  let fixture: ComponentFixture<Preprocess>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [Preprocess]
    })
    .compileComponents();

    fixture = TestBed.createComponent(Preprocess);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
