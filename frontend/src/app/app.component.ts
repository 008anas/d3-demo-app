import { Component, OnInit } from '@angular/core';
import { Meta } from '@angular/platform-browser';

import { TitleService } from '@services/title.service';
import { main } from '@config/main';

@Component({
  selector: 'sqy-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.scss'],
  providers: [TitleService]
})
export class AppComponent implements OnInit {

  constructor(
    private metaTagSrvc: Meta,
    private titleSrvc: TitleService
  ) { }

  ngOnInit() {
    this.metaTagSrvc.addTags([
      { name: 'keywords', content: 'Sequence Optimizator' },
      { name: 'keywords', content: main.appName },
      { name: 'robots', content: 'index, follow' },
      { name: 'author', content: 'Anas Gharrab' },
      { name: 'viewport', content: 'width=device-width, initial-scale=1' },
      { name: 'date', content: new Date().toString() },
      { charset: 'UTF-8' }
    ]);
    this.titleSrvc.init();
  }

}
