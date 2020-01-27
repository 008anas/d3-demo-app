import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { WorkspaceComponent } from './workspace.component';
import { HistoryComponent } from './history/history.component';
import { HistoryResolver } from './shared/history.resolver';
import { routes as rts } from '../config/routes';

const routes: Routes = [
  { path: '', component: WorkspaceComponent },
  { path: rts.workspace.detail, component: HistoryComponent, resolve: { history: HistoryResolver } }
];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule],
  providers: [HistoryResolver]
})
export class WorkspaceRoutingModule { }
