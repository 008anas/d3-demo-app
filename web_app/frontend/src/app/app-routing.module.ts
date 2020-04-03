import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { routes as rts } from '@config/routes';
import { HomeComponent } from './home/home.component';
import { ContactComponent } from './contact/contact.component';
import { Error404Component } from './core/errors/error404/error404.component';
import { Error500Component } from './core/errors/error500/error500.component';
import { ConstructResolver } from '@resolvers/construct.resolver';

const routes: Routes = [
  { path: '', redirectTo: '/', pathMatch: 'full' },
  { path: '', component: HomeComponent, data: { title: 'Home' } },
  { path: rts.home, component: HomeComponent, data: { title: 'Home' } },
  { path: rts.workspace.root, loadChildren: () => import('./workspace/workspace.module').then(m => m.WorkspaceModule) },
  { path: rts.optimize.root, loadChildren: () => import('./optimizer/optimizer.module').then(m => m.OptimizerModule) },
  { path: rts.vector.root, loadChildren: () => import('./vector-editor/vector-editor.module').then(m => m.VectorEditorModule) },
  { path: rts.contact, component: ContactComponent, data: { title: 'Contact us' } },
  { path: rts.error404, component: Error404Component, data: { title: 'Page not found' } },
  { path: rts.error500, component: Error500Component, data: { title: 'Internal Server Error' } },

  // Otherwise redirect to 404
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
