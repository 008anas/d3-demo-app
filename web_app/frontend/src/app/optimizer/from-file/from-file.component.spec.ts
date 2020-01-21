import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { FromFileComponent } from './from-file.component';

describe('FromFileComponent', () => {
  let component: FromFileComponent;
  let fixture: ComponentFixture<FromFileComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ FromFileComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(FromFileComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
