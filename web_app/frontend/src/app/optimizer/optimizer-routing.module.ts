import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { routes as rts } from '@config/routes';
import { SketcherComponent } from './sketcher/sketcher.component';
import { FromFileComponent } from './from-file/from-file.component';
import { CanDeactivateGuard } from '@services/can-deactivate-guard.service';

const routes: Routes = [
  { path: '', component: SketcherComponent, data: { title: 'Construct sketcher' } },
  { path: rts.optimize.sketcher, component: SketcherComponent, data: { title: 'Construct sketcher' }, canDeactivate: [CanDeactivateGuard] },
  { path: rts.optimize.from_file, component: FromFileComponent, data: { title: 'Upload file' } }
];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule],
  providers: [CanDeactivateGuard]
})
export class OptimizerRoutingModule { }
