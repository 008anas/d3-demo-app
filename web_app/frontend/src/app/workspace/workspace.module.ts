import { NgModule, CUSTOM_ELEMENTS_SCHEMA } from '@angular/core';
import { CommonModule } from '@angular/common';

import { WorkspaceRoutingModule } from './workspace-routing.module';
import { WorkspaceComponent } from './workspace.component';
import { HistoryComponent } from './history/history.component';
import { HistoryLoadingComponent } from './shared/history-loading/history-loading.component';


@NgModule({
  declarations: [
    WorkspaceComponent,
    HistoryComponent,
    HistoryLoadingComponent
  ],
  imports: [
    CommonModule,
    WorkspaceRoutingModule
  ],
  schemas: [
    CUSTOM_ELEMENTS_SCHEMA
  ]
})
export class WorkspaceModule { }
