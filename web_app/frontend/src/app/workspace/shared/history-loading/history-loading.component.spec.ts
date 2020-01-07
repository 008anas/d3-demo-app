import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { HistoryLoadingComponent } from './history-loading.component';

describe('HistoryLoadingComponent', () => {
  let component: HistoryLoadingComponent;
  let fixture: ComponentFixture<HistoryLoadingComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ HistoryLoadingComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(HistoryLoadingComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
