import { NgModule } from '@angular/core';
import { CommonModule } from '@angular/common';

import { NzModalModule } from 'ng-zorro-antd/modal';
import { NzSliderModule } from 'ng-zorro-antd/slider';
import { NzDropDownModule } from 'ng-zorro-antd/dropdown';

import { OptimizerRoutingModule } from './optimizer-routing.module';
import { OptimizerComponent } from './optimizer.component';
import { SketcherComponent } from './sketcher/sketcher.component';
import { SharedModule } from '../shared/shared.module';
import { TrackDetailsComponent } from './track-details/track-details.component';
import { FilterTracksPipe } from './shared/filter-tracks.pipe';


@NgModule({
  declarations: [
    OptimizerComponent,
    SketcherComponent,
    TrackDetailsComponent,
    FilterTracksPipe
  ],
  imports: [
    CommonModule,
    OptimizerRoutingModule,
    SharedModule,
    NzModalModule,
    NzSliderModule,
    NzDropDownModule
  ]
})
export class OptimizerModule { }
