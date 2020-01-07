import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClientModule } from '@angular/common/http';

import { ErrorsListComponent } from './errors-list/errors-list.component';
import { SidebarComponent } from './sidebar/sidebar.component';
import { ColorPickerComponent } from './color-picker/color-picker.component';

@NgModule({
  declarations: [
    ErrorsListComponent,
    SidebarComponent,
    ColorPickerComponent
  ],
  imports: [
    CommonModule,
    FormsModule,
    HttpClientModule
  ],
  exports: [
    ErrorsListComponent,
    SidebarComponent,
    ColorPickerComponent,
    FormsModule
  ]
})
export class SharedModule { }
