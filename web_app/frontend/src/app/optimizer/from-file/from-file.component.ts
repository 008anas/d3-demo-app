import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';
import { ActivatedRoute, Router } from '@angular/router';
import { HttpClient } from '@angular/common/http';
import { finalize } from 'rxjs/operators';

import { environment as env } from '@env/environment';

import { Specie } from '@models/specie';
import { Construct } from '@models/construct';
import { SpecieService } from '@services/specie.service';
import { UserHistory } from 'app/workspace/shared/user-history';

@Component({
  selector: 'sqy-from-file',
  templateUrl: './from-file.component.html',
  styleUrls: ['./from-file.component.scss']
})
export class FromFileComponent implements OnInit, OnDestroy {

  response: any = null;
  sub: Subscription;
  specie: Specie = new Specie();
  species: Specie[] = [];
  isLoading = false;
  history: UserHistory = null;
  isUploading = false;
  construct: Construct = null;
  exampleFile = 'assets/files/NC_044937.gbk';
  endpoint: string;

  constructor(
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private router: Router,
    private http: HttpClient
  ) {
    this.endpoint = env.endpoints.api + '/constructs/from-genbank';
  }

  ngOnInit() {
    this.sub = this.route.queryParams.subscribe(params => this.specie.slug = params.specie || null);
    this.getSpecie();
    this.getSpecies();
  }

  ngOnDestroy() {
    if (this.sub) { this.sub.unsubscribe(); }
  }

  getSpecie() {
    this.isLoading = true;
    this.specieSrvc.getBySlug(this.specie)
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(data => this.specie.deserialize(data));
  }

  getSpecies() {
    this.isLoading = true;
    this.specieSrvc.getAll()
      .pipe(finalize(() => this.isLoading = false))
      .subscribe(
        data =>
          this.species = data.map((e: any) => {
            const specie = new Specie().deserialize(e);
            if (!this.specie.slug && specie.default) {
              this.specie = Object.assign({}, specie);
            }
            return specie;
          })
      );
  }

  handleUploaded(construct: Construct) {
    this.construct = new Construct().deserialize(construct);
  }

  submit() {
    this.response = null;
    this.construct["specie_tax_id"] = this.specie.tax_id;
    this.http
      .post(`${env.endpoints.api}/optimize_seq/from-sketch`, this.construct)
      .subscribe(
        (data: UserHistory) => {
          this.history = new UserHistory().deserialize(data);
          setTimeout(() => {
            this.router.navigate(["/workspace", this.history.id]);
          }, 3000);
        },
        err => {
          // this.notify.error(err.msg, "bottom-right");
        }
      );
  }

}
