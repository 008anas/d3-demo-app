import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { routes as rts } from '@config/routes';
import { HomeComponent } from './home/home.component';
import { DocumentationComponent } from './documentation/documentation.component';
import { ContactComponent } from './contact/contact.component';
import { Error404Component } from './core/errors/error404/error404.component';
import { Error500Component } from './core/errors/error500/error500.component';
import { ConstructResolver } from './construct/shared/construct.resolver';
import { VectorEditorComponent } from './vector-editor/vector-editor.component';

const routes: Routes = [
  { path: '', redirectTo: '/', pathMatch: 'full' },
  { path: '', component: HomeComponent, data: { title: 'Select specie' } },
  { path: rts.workspace.root, loadChildren: () => import('./workspace/workspace.module').then(m => m.WorkspaceModule) },
  { path: rts.optimize.root, loadChildren: () => import('./optimizer/optimizer.module').then(m => m.OptimizerModule) },
  // { path: rts.construct.root, loadChildren: () => import('./construct/construct.module').then(m => m.ConstructModule) },
  { path: 'editor', component: VectorEditorComponent, data: { title: 'Editor' } },
  { path: 'editor/:construct', component: VectorEditorComponent, resolve: { construct: ConstructResolver }, data: { title: 'Editor' } },
  { path: rts.documentation, component: DocumentationComponent, data: { title: 'Documentation' } },
  { path: rts.contact, component: ContactComponent, data: { title: 'Contact us' } },
  { path: rts.error404, component: Error404Component, data: { title: 'Request page not found 404' } },
  { path: rts.error500, component: Error500Component, data: { title: 'Internal Server Error 500' } },

  // otherwise redirect to 404
  { path: '**', redirectTo: '/' + rts.error404 }
];

@NgModule({
  imports: [RouterModule.forRoot(routes, {
    onSameUrlNavigation: 'reload'
  })],
  exports: [RouterModule],
  providers: [ConstructResolver]
})
export class AppRoutingModule { }
