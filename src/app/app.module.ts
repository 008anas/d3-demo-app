import { BrowserModule } from '@angular/platform-browser';
import { NgModule, ErrorHandler } from '@angular/core';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { HTTP_INTERCEPTORS, HttpClientModule } from '@angular/common/http';

import { NgProgressModule } from '@ngx-progressbar/core';
import { NgProgressHttpModule } from '@ngx-progressbar/http';
import { NZ_I18N, en_US } from 'ng-zorro-antd';
import { NzSelectModule } from 'ng-zorro-antd/select';
import { NzMessageModule } from 'ng-zorro-antd/message';

import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { CoreModule } from './core/core.module';
import { HomeComponent } from './home/home.component';
import { ContactComponent } from './contact/contact.component';
import { SharedModule } from './shared/shared.module';
import { ErrorInterceptor } from './core/errors/error.interceptor';
import { AppHttpInterceptor } from './core/app.interceptor';
import { SentryErrorHandler } from './core/errors/sentry-error-handler.service';
import { LoaderComponent } from '@components/loader/loader.component';
import { VectorEditorModule } from './vector-editor/vector-editor.module';
import { SuccessComponent } from './contact/success/success.component';

@NgModule({
  declarations: [
    AppComponent,
    HomeComponent,
    ContactComponent,
    LoaderComponent,
    SuccessComponent
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
    NzMessageModule,
    VectorEditorModule
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
