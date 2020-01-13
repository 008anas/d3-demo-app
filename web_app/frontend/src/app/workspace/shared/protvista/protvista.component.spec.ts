import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { ProtvistaComponent } from './protvista.component';

describe('ProtvistaComponent', () => {
  let component: ProtvistaComponent;
  let fixture: ComponentFixture<ProtvistaComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ ProtvistaComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(ProtvistaComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
