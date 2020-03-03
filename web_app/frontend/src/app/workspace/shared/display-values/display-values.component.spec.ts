import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { DisplayValuesComponent } from './display-values.component';

describe('DisplayValuesComponent', () => {
  let component: DisplayValuesComponent;
  let fixture: ComponentFixture<DisplayValuesComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ DisplayValuesComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(DisplayValuesComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
