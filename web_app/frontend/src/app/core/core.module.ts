import { NgModule, Optional, SkipSelf } from '@angular/core';
import { CommonModule } from '@angular/common';
import { RouterModule } from '@angular/router';

import { CookieService } from 'ngx-cookie-service';
import { NzDropDownModule } from 'ng-zorro-antd/dropdown';

import { throwIfAlreadyLoaded } from './module-import-guard';
import { FooterComponent } from './footer/footer.component';
import { NavComponent } from './nav/nav.component';
import { Error404Component } from './errors/error404/error404.component';
import { Error500Component } from './errors/error500/error500.component';


@NgModule({
  declarations: [
    FooterComponent,
    NavComponent,
    Error404Component,
    Error500Component
  ],
  imports: [
    CommonModule,
    RouterModule,
    NzDropDownModule
  ],
  exports: [
    FooterComponent,
    NavComponent
  ],
  providers: [CookieService]
})
export class CoreModule {
  constructor(@Optional() @SkipSelf() parentModule: CoreModule) {
    throwIfAlreadyLoaded(parentModule, 'CoreModule');
  }
}
