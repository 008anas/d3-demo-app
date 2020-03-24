import { Component, OnInit, OnDestroy } from '@angular/core';
import { Router, ActivatedRoute } from '@angular/router';
import { Subscription } from 'rxjs/internal/Subscription';
import { finalize } from 'rxjs/operators';

import { Specie } from '@models/specie';
import { SpecieService } from '@services/specie.service';
import { HistoryService } from 'app/workspace/shared/history.service';
import { UserHistory } from 'app/workspace/shared/user-history';
import Utils from 'app/shared/utils';

@Component({
  selector: 'sqy-home',
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.scss']
})
export class HomeComponent implements OnInit, OnDestroy {

  sub: Subscription;
  specie: Specie = null;
  species: Specie[] = [];
  history_id: string = null;
  isLoading = false;
  response: any = null;
  isRetriving = false;
  option = 1;
  uuidRegex = Utils.uuidRegex;
  carousel = [
    'https://www.earlham.ac.uk/sites/default/files/images/articles/Synthetic%20Biology%20Spotlight%20Yaomin%20Cai/synthetic-biology-spotlight-yaomin-cai-hero.jpg',
    'https://www.genengnews.com/wp-content/uploads/2019/09/Sep3_2019_Getty_157647511_DNASequence-1-1068x712.jpg',
    'https://sitechecker.pro/wp-content/uploads/2017/12/search-engines.png'
  ];

  constructor(
    private specieSrvc: SpecieService,
    private route: ActivatedRoute,
    private router: Router,
    private historySrvc: HistoryService
  ) { }

  ngOnInit() {
    this.isLoading = true;
    this.sub = this.route.queryParams.subscribe(params => {
      if (params.option && params.option <= this.carousel.length) {
        this.option = Number.parseInt(params.option, 0);
      }
    });
    this.getSpecies();
  }

  ngOnDestroy() {
    if (this.sub) {
      this.sub.unsubscribe();
    }
  }

  getSpecies() {
    this.isLoading = true;
    this.specieSrvc.getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(
        data => {
          this.species = data.map((e: any) => {
            const specie = new Specie().deserialize(e);
            if (specie.default) { this.specie = specie; }
            return specie;
          });
        }
      );
  }

  getHistory() {
    if (this.history_id) {
      this.isRetriving = true;
      this.response = null;
      this.historySrvc.getByIdNot404(this.history_id)
        .pipe(finalize(() => this.isRetriving = false))
        .subscribe(
          (data: UserHistory) => {
            this.router.navigate(['/workspace', data.id]);
          },
          err => this.response = err
        );
    }
  }

  goTo(path: string) {
    if (this.specie) {
      this.router.navigate(['/optimize', path], { queryParams: { specie: this.specie.slug } });
    } else {
      alert('First choose a specie to proceed');
    }
  }
}
