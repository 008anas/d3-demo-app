import { Component, OnInit } from '@angular/core';
import { Router } from '@angular/router';
import { finalize } from 'rxjs/operators';

import { Specie } from '../shared/models/specie';
import { SpecieService } from '../shared/services/specie.service';

@Component({
  selector: 'sqy-home',
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.scss']
})
export class HomeComponent implements OnInit {

  specie: Specie = null;
  species: Specie[] = [];
  history_id: string = null;
  isLoading = false;

  constructor(
    private specieSrvc: SpecieService,
    private router: Router
  ) { }

  ngOnInit() {
    this.getSpecies();
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

  goTo(path: string) {
    this.specie ? this.router.navigate(['/optimize', path], { queryParams: { specie: this.specie.slug } }) : alert('First choose a specie to proceed');
  }
}
