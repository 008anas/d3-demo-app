import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';

import { ConstructRoutingModule } from './construct-routing.module';
import { ConstructComponent } from './construct.component';


@NgModule({
  declarations: [ConstructComponent],
  imports: [
    CommonModule,
    ConstructRoutingModule
  ]
})
export class ConstructModule { }
