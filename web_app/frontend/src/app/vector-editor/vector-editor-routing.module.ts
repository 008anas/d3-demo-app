import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { routes as rts } from '@config/routes';
import { VectorComponent } from './vector/vector.component';
import { HistoryResolver } from 'app/workspace/shared/history.resolver';

const routes: Routes = [
  { path: '', component: VectorComponent, data: { title: 'Vector editor' } },
  { path: rts.editor.withHistory, component: VectorComponent, resolve: { history: HistoryResolver } }
];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule],
  providers: [HistoryResolver]
})
export class VectorEditorRoutingModule { }
