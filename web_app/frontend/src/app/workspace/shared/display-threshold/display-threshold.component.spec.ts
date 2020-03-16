import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { DisplayThresholdComponent } from './display-threshold.component';

describe('DisplayThresholdComponent', () => {
  let component: DisplayThresholdComponent;
  let fixture: ComponentFixture<DisplayThresholdComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ DisplayThresholdComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(DisplayThresholdComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
