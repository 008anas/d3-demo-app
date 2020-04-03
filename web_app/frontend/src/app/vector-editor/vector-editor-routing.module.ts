import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { routes as rts } from '@config/routes';
import { HistoryResolver } from 'app/workspace/shared/history.resolver';
import { EditorComponent } from './editor/editor.component';

const routes: Routes = [
  { path: '', component: EditorComponent, data: { title: 'Vector editor' } },
  { path: rts.vector.withHistory, component: EditorComponent, resolve: { history: HistoryResolver } }
];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule],
  providers: [HistoryResolver]
})
export class VectorEditorRoutingModule { }
