import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';

import { VectorEditorRoutingModule } from './vector-editor-routing.module';
import { VectorComponent } from './vector/vector.component';


@NgModule({
  declarations: [
    VectorComponent
  ],
  imports: [
    CommonModule,
    VectorEditorRoutingModule
  ]
})
export class VectorEditorModule { }
