import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { TracksPickerComponent } from './tracks-picker.component';

describe('TracksPickerComponent', () => {
  let component: TracksPickerComponent;
  let fixture: ComponentFixture<TracksPickerComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ TracksPickerComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(TracksPickerComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
