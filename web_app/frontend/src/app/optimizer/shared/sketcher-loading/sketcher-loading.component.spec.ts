import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { SketcherLoadingComponent } from './sketcher-loading.component';

describe('SketcherLoadingComponent', () => {
  let component: SketcherLoadingComponent;
  let fixture: ComponentFixture<SketcherLoadingComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ SketcherLoadingComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SketcherLoadingComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
