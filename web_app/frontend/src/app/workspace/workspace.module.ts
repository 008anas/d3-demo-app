import { NgModule, CUSTOM_ELEMENTS_SCHEMA } from '@angular/core';
import { CommonModule } from '@angular/common';
import { ReactiveFormsModule } from '@angular/forms';

import { NzModalModule } from 'ng-zorro-antd/modal';
import { NzSelectModule } from 'ng-zorro-antd/select';
import { NzDropDownModule } from 'ng-zorro-antd/dropdown';
import { NzPopoverModule } from 'ng-zorro-antd/popover';

import { WorkspaceRoutingModule } from './workspace-routing.module';
import { WorkspaceComponent } from './workspace.component';
import { HistoryComponent } from './history/history.component';
import { HistoryLoadingComponent } from './shared/history-loading/history-loading.component';
import { SharedModule } from '../shared/shared.module';
import { ResultsViewerComponent } from './shared/results-viewer/results-viewer.component';
import { DisplayValuesComponent } from './shared/display-values/display-values.component';
import { DisplayThresholdComponent } from './shared/display-threshold/display-threshold.component';


@NgModule({
  declarations: [
    WorkspaceComponent,
    HistoryComponent,
    HistoryLoadingComponent,
    ResultsViewerComponent,
    DisplayValuesComponent,
    DisplayThresholdComponent
  ],
  imports: [
    CommonModule,
    WorkspaceRoutingModule,
    SharedModule,
    NzModalModule,
    NzSelectModule,
    NzDropDownModule,
    NzPopoverModule,
    ReactiveFormsModule
  ],
  schemas: [
    CUSTOM_ELEMENTS_SCHEMA
  ]
})
export class WorkspaceModule { }
