import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';

import { NzModalModule } from 'ng-zorro-antd/modal';
import { NzSliderModule } from 'ng-zorro-antd/slider';
import { NzDropDownModule } from 'ng-zorro-antd/dropdown';
import { NzUploadModule } from 'ng-zorro-antd/upload';
import { NzAlertModule } from 'ng-zorro-antd/alert';
import { NzDrawerModule } from 'ng-zorro-antd/drawer';

import { OptimizerRoutingModule } from './optimizer-routing.module';
import { OptimizerComponent } from './optimizer.component';
import { SketcherComponent } from './sketcher/sketcher.component';
import { SharedModule } from '../shared/shared.module';
import { TrackDetailsComponent } from './track-details/track-details.component';
import { FilterTracksPipe } from './shared/filter-tracks.pipe';
import { FromFileComponent } from './from-file/from-file.component';
import { ReactiveFormsModule } from '@angular/forms';
import { TracksPickerComponent } from './shared/tracks-picker/tracks-picker.component';


@NgModule({
  declarations: [
    OptimizerComponent,
    SketcherComponent,
    TrackDetailsComponent,
    FilterTracksPipe,
    FromFileComponent,
    TracksPickerComponent
  ],
  imports: [
    CommonModule,
    OptimizerRoutingModule,
    SharedModule,
    NzModalModule,
    NzSliderModule,
    NzDropDownModule,
    ReactiveFormsModule,
    NzUploadModule,
    NzAlertModule,
    NzDrawerModule
  ]
})
export class OptimizerModule { }
