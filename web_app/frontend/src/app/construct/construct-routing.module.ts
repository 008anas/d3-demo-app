import { NgModule } from '@angular/core';
import { Routes, RouterModule } from '@angular/router';

import { ConstructComponent } from './construct.component';

const routes: Routes = [{ path: '', component: ConstructComponent }];

@NgModule({
  imports: [RouterModule.forChild(routes)],
  exports: [RouterModule]
})
export class ConstructRoutingModule { }
