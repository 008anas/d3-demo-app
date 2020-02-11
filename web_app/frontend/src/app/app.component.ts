import { Component, Inject, PLATFORM_ID, OnInit } from '@angular/core';
import { Meta } from '@angular/platform-browser';
import { isPlatformBrowser } from '@angular/common';

import { TitleService } from './shared/services/title.service';

@Component({
  selector: 'sqy-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.scss'],
  providers: [TitleService]
})
export class AppComponent implements OnInit {

  isOnline: boolean;

  constructor(
    private metaTagSrvc: Meta,
    private titleSrvc: TitleService,
    @Inject(PLATFORM_ID) private platformId: object
  ) {
    isPlatformBrowser(this.platformId) ? this.isOnline = navigator.onLine : this.isOnline = true;
  }

  ngOnInit() {
    this.metaTagSrvc.addTags([
      { name: 'keywords', content: 'Sequence Optimizator' },
      { name: 'keywords', content: 'SQRUTINY' },
      { name: 'robots', content: 'index, follow' },
      { name: 'author', content: 'Anas Gharrab' },
      { name: 'viewport', content: 'width=device-width, initial-scale=1' },
      { name: 'date', content: '2019-10-31', scheme: 'YYYY-MM-DD' },
      { charset: 'UTF-8' }
    ]);
    this.titleSrvc.init();
  }

}
