import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { environment as env } from 'src/environments/environment';
import { SketcherComponent } from './sketcher/sketcher.component';
import { FromFileComponent } from './from-file/from-file.component';

const routes: Routes = [
  { path: '', component: SketcherComponent, data: { title: 'Construct sketcher' } },
  { path: env.routes.optimize.sketcher, component: SketcherComponent, data: { title: 'Construct sketcher' } },
  { path: env.routes.optimize.from_file, component: FromFileComponent, data: { title: 'Upload file' } }
];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule]
})
export class OptimizerRoutingModule { }
