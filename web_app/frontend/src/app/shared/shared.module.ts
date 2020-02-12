import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

import { SafeDatePipe } from './pipes/date.pipe';
import { CopyClipboardDirective } from './directives/copy-clipboard.directive';
import { ErrorsListComponent } from './components/errors-list/errors-list.component';
import { SidebarComponent } from './components/sidebar/sidebar.component';
import { ColorPickerComponent } from './components/color-picker/color-picker.component';
import { IfNotIconComponent } from './components/if-not-icon/if-not-icon.component';
import { ProgressBarComponent } from './components/progress-bar/progress-bar.component';
import { FileUploaderComponent } from './components/file-uploader/file-uploader.component';

@NgModule({
  declarations: [
    ErrorsListComponent,
    SidebarComponent,
    ColorPickerComponent,
    SafeDatePipe,
    CopyClipboardDirective,
    IfNotIconComponent,
    ProgressBarComponent,
    FileUploaderComponent
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
    CopyClipboardDirective,
    IfNotIconComponent,
    ProgressBarComponent
  ]
})
export class SharedModule { }
