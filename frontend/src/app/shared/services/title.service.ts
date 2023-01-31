import { Injectable } from '@angular/core';

import { Router, NavigationEnd, ActivatedRoute } from '@angular/router';
import { Title } from '@angular/platform-browser';
import { filter, map } from 'rxjs/operators';
import { main } from '@config/main';

const SEPARATOR = ' | ';

@Injectable()
export class TitleService {

  static ucFirst(value: string) {
    if (!value) { return value; }
    return value.charAt(0).toUpperCase() + value.slice(1);
  }

  constructor(
    private router: Router,
    private activatedRoute: ActivatedRoute,
    private titleSrvc: Title
  ) { }

  init() {
    this.router.events.pipe(
      filter((event) => event instanceof NavigationEnd),
      map(() => {
        let route = this.activatedRoute;
        while (route.firstChild) { route = route.firstChild; }
        return route;
      }),
      filter((route) => route.outlet === 'primary'),
      map((route) => route.snapshot),
      map((snapshot) => {
        if (snapshot.data.title) {
          if (snapshot.paramMap.get('id') !== null) {
            return snapshot.data.title + SEPARATOR + snapshot.paramMap.get('id');
          }
          return snapshot.data.title;
        } else {
          // If not, we do a little magic on the url to create an approximation
          return this.router.url.split('/').reduce((acc, frag) => {
            if (acc && frag) { acc += SEPARATOR; }
            return acc + TitleService.ucFirst(frag);
          });
        }
      })
    ).subscribe((pathStr) => this.titleSrvc.setTitle(`${pathStr}${SEPARATOR}${main.appName}`));
  }

  setTitle(str: string) {
    this.titleSrvc.setTitle(`${str}${SEPARATOR}${main.appName}`);
  }
}
