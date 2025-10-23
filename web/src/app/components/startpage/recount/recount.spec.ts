import { ComponentFixture, TestBed } from '@angular/core/testing';

import { Recount } from './recount';

describe('Recount', () => {
  let component: Recount;
  let fixture: ComponentFixture<Recount>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [Recount]
    })
    .compileComponents();

    fixture = TestBed.createComponent(Recount);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
