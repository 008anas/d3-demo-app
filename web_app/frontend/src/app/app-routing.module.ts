import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { environment as env } from 'src/environments/environment';
import { HomeComponent } from './home/home.component';
import { DocumentationComponent } from './documentation/documentation.component';
import { ContactComponent } from './contact/contact.component';
import { Error404Component } from './core/errors/error404/error404.component';
import { Error500Component } from './core/errors/error500/error500.component';
import { EditorComponent } from './editor/editor.component';
import { ConstructResolver } from './construct/shared/construct.resolver';

const routes: Routes = [
  { path: '', redirectTo: '/', pathMatch: 'full' },
  { path: '', component: HomeComponent, data: { title: 'Select specie' } },
  { path: env.routes.workspace.root, loadChildren: () => import('./workspace/workspace.module').then(m => m.WorkspaceModule) },
  { path: env.routes.optimize.root, loadChildren: () => import('./optimizer/optimizer.module').then(m => m.OptimizerModule) },
  { path: env.routes.construct.root, loadChildren: () => import('./construct/construct.module').then(m => m.ConstructModule) },
  { path: 'editor', component: EditorComponent, data: { title: 'Editor' } },
  { path: 'editor/:construct', component: EditorComponent, resolve: { construct: ConstructResolver }, data: { title: 'Editor' } },
  { path: env.routes.documentation, component: DocumentationComponent, data: { title: 'Documentation' } },
  { path: env.routes.contact, component: ContactComponent, data: { title: 'Contact us' } },
  { path: env.routes.error404, component: Error404Component, data: { title: 'Request page not found 404' } },
  { path: env.routes.error500, component: Error500Component, data: { title: 'Internal Server Error 500' } },

  // otherwise redirect to 404
  { path: '**', redirectTo: '/' + env.routes.error404 }
];

@NgModule({
  imports: [RouterModule.forRoot(routes, {
    onSameUrlNavigation: 'reload'
  })],
  exports: [RouterModule],
  providers: [ConstructResolver]
})
export class AppRoutingModule { }
