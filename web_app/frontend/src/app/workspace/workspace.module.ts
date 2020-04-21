import { NgModule, CUSTOM_ELEMENTS_SCHEMA } from '@angular/core';
import { CommonModule } from '@angular/common';
import { ReactiveFormsModule } from '@angular/forms';

import { NzModalModule } from 'ng-zorro-antd/modal';
import { NzSelectModule } from 'ng-zorro-antd/select';
import { NzDropDownModule } from 'ng-zorro-antd/dropdown';

import { WorkspaceRoutingModule } from './workspace-routing.module';
import { WorkspaceComponent } from './workspace.component';
import { HistoryComponent } from './history/history.component';
import { HistoryLoadingComponent } from './shared/history-loading/history-loading.component';
import { SharedModule } from '../shared/shared.module';
import { ResultViewerComponent } from './shared/result-viewer/result-viewer.component';
import { DisplayValuesComponent } from './shared/display-values/display-values.component';
import { SetHighlightComponent } from './shared/set-highlight/set-highlight.component';
import { ExportModalComponent } from './shared/export-modal/export-modal.component';
import { SetCutoffComponent } from './shared/set-cutoff/set-cutoff.component';

@NgModule({
  declarations: [
    WorkspaceComponent,
    HistoryComponent,
    HistoryLoadingComponent,
    ResultViewerComponent,
    DisplayValuesComponent,
    SetHighlightComponent,
    ExportModalComponent,
    SetCutoffComponent
  ],
  imports: [
    CommonModule,
    WorkspaceRoutingModule,
    SharedModule,
    NzModalModule,
    NzSelectModule,
    NzDropDownModule,
    ReactiveFormsModule
  ],
  schemas: [
    CUSTOM_ELEMENTS_SCHEMA
  ]
})
export class WorkspaceModule { }
