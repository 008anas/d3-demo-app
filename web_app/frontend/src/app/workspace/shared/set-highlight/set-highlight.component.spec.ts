import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { SetHighlightComponent } from './set-highlight.component';

describe('SetHighlightComponent', () => {
  let component: SetHighlightComponent;
  let fixture: ComponentFixture<SetHighlightComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ SetHighlightComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SetHighlightComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
