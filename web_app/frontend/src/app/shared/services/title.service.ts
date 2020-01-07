import { Injectable } from '@angular/core';

import { Router, NavigationEnd, ActivatedRoute } from '@angular/router';
import { Title } from '@angular/platform-browser';
import { filter, map } from 'rxjs/operators';
import { environment as env } from 'src/environments/environment';

const SEPARATOR = ' | ';

@Injectable()
export class TitleService {
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
    ).subscribe((pathStr) => this.titleSrvc.setTitle(`${pathStr}${SEPARATOR}${env.name}`));
  }

  setTitle(str: string) {
    this.titleSrvc.setTitle(`${str}${SEPARATOR}${env.name}`);
  }

  static ucFirst(string: string) {
    if (!string) { return string; }
    return string.charAt(0).toUpperCase() + string.slice(1);
  }
}
