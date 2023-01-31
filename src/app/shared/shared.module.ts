import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

import { NzModalModule } from 'ng-zorro-antd/modal';
import { NzMessageModule } from 'ng-zorro-antd/message';
import { NzPopoverModule } from 'ng-zorro-antd/popover';

import { SafeDatePipe } from '@pipes/date.pipe';
import { CopyClipboardDirective } from '@directives/copy-clipboard.directive';
import { ErrorsListComponent } from '@components/errors-list/errors-list.component';
import { IfNotIconComponent } from '@components/if-not-icon/if-not-icon.component';
import { TextModalComponent } from '@components/text-modal/text-modal.component';
import { InfoPopupComponent } from './components/info-popup/info-popup.component';
import { TrackDisplayComponent } from './components/track-display/track-display.component';

@NgModule({
  declarations: [
    ErrorsListComponent,
    SafeDatePipe,
    CopyClipboardDirective,
    IfNotIconComponent,
    TextModalComponent,
    InfoPopupComponent,
    TrackDisplayComponent
  ],
  imports: [
    CommonModule,
    FormsModule,
    NzModalModule,
    NzMessageModule,
    NzPopoverModule
  ],
  exports: [
    FormsModule,
    ErrorsListComponent,
    SafeDatePipe,
    CopyClipboardDirective,
    IfNotIconComponent,
    TextModalComponent,
    InfoPopupComponent,
    NzPopoverModule,
    TrackDisplayComponent
  ]
})
export class SharedModule { }
