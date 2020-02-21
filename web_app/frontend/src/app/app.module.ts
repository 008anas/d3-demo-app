import { BrowserModule } from '@angular/platform-browser';
import { NgModule, ErrorHandler } from '@angular/core';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { HTTP_INTERCEPTORS, HttpClientModule } from '@angular/common/http';

import { NgProgressModule } from '@ngx-progressbar/core';
import { NgProgressHttpModule } from '@ngx-progressbar/http';
import { NZ_I18N, en_US } from 'ng-zorro-antd';
import { NzSelectModule } from 'ng-zorro-antd/select';
import { NzAlertModule } from 'ng-zorro-antd/alert';
import { NzCarouselModule } from 'ng-zorro-antd/carousel';

import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { CoreModule } from './core/core.module';
import { HomeComponent } from './home/home.component';
import { ContactComponent } from './contact/contact.component';
import { SharedModule } from './shared/shared.module';
import { ErrorInterceptor } from './core/errors/error.interceptor';
import { AppHttpInterceptor } from './core/app.interceptor';
import { SentryErrorHandler } from './core/errors/sentry-error-handler.service';
import { VectorEditorComponent } from './vector-editor/vector-editor.component';
import { NotifyComponent } from './shared/components/notify/notify.component';

@NgModule({
  declarations: [
    AppComponent,
    HomeComponent,
    NotifyComponent,
    ContactComponent,
    VectorEditorComponent
  ],
  imports: [
    BrowserModule,
    BrowserAnimationsModule,
    AppRoutingModule,
    NgProgressHttpModule,
    NgProgressModule.withConfig({
      color: '#fff'
    }),
    HttpClientModule,
    CoreModule,
    SharedModule,
    NzSelectModule,
    NzAlertModule,
    NzCarouselModule
  ],
  providers: [
    { provide: NZ_I18N, useValue: en_US },
    { provide: HTTP_INTERCEPTORS, useClass: ErrorInterceptor, multi: true },
    { provide: HTTP_INTERCEPTORS, useClass: AppHttpInterceptor, multi: true },
    { provide: ErrorHandler, useClass: SentryErrorHandler }
  ],
  bootstrap: [AppComponent]
})
export class AppModule { }
