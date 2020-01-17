import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { IfNotIconComponent } from './if-not-icon.component';

describe('IfNotIconComponent', () => {
  let component: IfNotIconComponent;
  let fixture: ComponentFixture<IfNotIconComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ IfNotIconComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(IfNotIconComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
