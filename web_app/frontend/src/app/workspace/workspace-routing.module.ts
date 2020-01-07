import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { environment as env } from 'src/environments/environment';

import { WorkspaceComponent } from './workspace.component';
import { HistoryComponent } from './history/history.component';
import { HistoryResolver } from './shared/history.resolver';

const routes: Routes = [
  { path: '', component: WorkspaceComponent },
  { path: env.routes.workspace.detail, component: HistoryComponent, resolve: { history: HistoryResolver } }
];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule],
  providers: [HistoryResolver]
})
export class WorkspaceRoutingModule { }
