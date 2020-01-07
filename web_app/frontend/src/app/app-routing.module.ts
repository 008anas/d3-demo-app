import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { environment as env } from 'src/environments/environment';
import { HomeComponent } from './home/home.component';
import { DocumentationComponent } from './documentation/documentation.component';
import { ContactComponent } from './contact/contact.component';
import { Error404Component } from './core/errors/error404/error404.component';
import { Error500Component } from './core/errors/error500/error500.component';


const routes: Routes = [
  { path: '', redirectTo: '/', pathMatch: 'full' },
  { path: 'workspace', loadChildren: () => import('./workspace/workspace.module').then(m => m.WorkspaceModule) },
  { path: '', component: HomeComponent, data: { title: 'Select specie' } },
  { path: env.routes.optimize.root, loadChildren: './optimizer/optimizer.module#OptimizerModule' },
  // // { path: env.routes.construct.root, loadChildren: './construct/construct.module#ConstructModule'},
  { path: env.routes.workspace.root, loadChildren: './workspace/workspace.module#WorkspaceModule' },
  { path: env.routes.documentation, component: DocumentationComponent, data: { title: 'Documentation' } },
  { path: env.routes.contact, component: ContactComponent, data: { title: 'Contact us' } },
  { path: env.routes.error404, component: Error404Component, data: { title: 'Request page not found 404' } },
  { path: env.routes.error500, component: Error500Component, data: { title: 'Internal Server Error 500' } },

  // otherwise redirect to 404
  { path: '**', redirectTo: '/' + env.routes.error404 },
  { path: 'optimizer', loadChildren: () => import('./optimizer/optimizer.module').then(m => m.OptimizerModule) },
  { path: 'construct', loadChildren: () => import('./construct/construct.module').then(m => m.ConstructModule) }
];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
