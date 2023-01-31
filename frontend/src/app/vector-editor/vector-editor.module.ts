import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';

import { VectorEditorRoutingModule } from './vector-editor-routing.module';
import { SharedModule } from 'app/shared/shared.module';
import { EditorComponent } from './editor/editor.component';
import { HistoryPickerComponent } from './shared/history-picker/history-picker.component';

@NgModule({
  declarations: [
    EditorComponent,
    HistoryPickerComponent
  ],
  imports: [
    CommonModule,
    VectorEditorRoutingModule,
    SharedModule
  ]
})
export class VectorEditorModule { }
