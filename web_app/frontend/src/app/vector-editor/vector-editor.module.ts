import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';

import { VectorEditorRoutingModule } from './vector-editor-routing.module';
import { VectorComponent } from './vector/vector.component';
import { SharedModule } from 'app/shared/shared.module';


@NgModule({
  declarations: [
    VectorComponent
  ],
  imports: [
    CommonModule,
    VectorEditorRoutingModule,
    SharedModule
  ]
})
export class VectorEditorModule { }
