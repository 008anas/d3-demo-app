import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { VectorEditorComponent } from './vector-editor.component';

describe('VectorEditorComponent', () => {
  let component: VectorEditorComponent;
  let fixture: ComponentFixture<VectorEditorComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ VectorEditorComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(VectorEditorComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
