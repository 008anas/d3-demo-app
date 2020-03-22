import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

import { NzModalModule } from 'ng-zorro-antd/modal';
import { NzMessageModule } from 'ng-zorro-antd/message';

import { SafeDatePipe } from '@pipes/date.pipe';
import { CopyClipboardDirective } from '@directives/copy-clipboard.directive';
import { ErrorsListComponent } from '@components/errors-list/errors-list.component';
import { IfNotIconComponent } from '@components/if-not-icon/if-not-icon.component';
import { BytesPipe } from '@pipes/bytes.pipe';
import { TextModalComponent } from '@components/text-modal/text-modal.component';

@NgModule({
  declarations: [
    ErrorsListComponent,
    SafeDatePipe,
    CopyClipboardDirective,
    IfNotIconComponent,
    BytesPipe,
    TextModalComponent
  ],
  imports: [
    CommonModule,
    FormsModule,
    NzModalModule,
    NzMessageModule
  ],
  exports: [
    FormsModule,
    ErrorsListComponent,
    SafeDatePipe,
    CopyClipboardDirective,
    IfNotIconComponent,
    TextModalComponent
  ]
})
export class SharedModule { }
