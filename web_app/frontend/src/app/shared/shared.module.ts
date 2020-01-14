import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

import { ErrorsListComponent } from './errors-list/errors-list.component';
import { SidebarComponent } from './sidebar/sidebar.component';
import { ColorPickerComponent } from './color-picker/color-picker.component';
import { SafeDatePipe } from './pipes/date.pipe';
import { CopyClipboardDirective } from './directives/copy-clipboard.directive';

@NgModule({
  declarations: [
    ErrorsListComponent,
    SidebarComponent,
    ColorPickerComponent,
    SafeDatePipe,
    CopyClipboardDirective
  ],
  imports: [
    CommonModule,
    FormsModule
  ],
  exports: [
    FormsModule,
    ErrorsListComponent,
    SidebarComponent,
    ColorPickerComponent,
    SafeDatePipe,
    CopyClipboardDirective
  ]
})
export class SharedModule { }
