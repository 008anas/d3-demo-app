import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { SetCutoffComponent } from './set-cutoff.component';

describe('SetCutoffComponent', () => {
  let component: SetCutoffComponent;
  let fixture: ComponentFixture<SetCutoffComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ SetCutoffComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SetCutoffComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
