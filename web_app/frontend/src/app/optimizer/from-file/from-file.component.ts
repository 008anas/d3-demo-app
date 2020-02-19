import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';
import { ActivatedRoute, Router } from '@angular/router';
import { HttpClient, HttpResponse } from '@angular/common/http';
import { finalize } from 'rxjs/operators';

import { UploadChangeParam } from 'ng-zorro-antd/upload';

import { environment as env } from '@env/environment';
import { Specie } from '@models/specie';
import { Construct } from '@models/construct';
import { SpecieService } from '@services/specie.service';
import { UserHistory } from 'app/workspace/shared/user-history';
import { NotifyService } from '@services/notify.service';

@Component({
  selector: 'sqy-from-file',
  templateUrl: './from-file.component.html',
  styleUrls: ['./from-file.component.scss']
})
export class FromFileComponent implements OnInit, OnDestroy {

  response: any = null;
  sub: Subscription;
  specie: Specie = null;
  species: Specie[] = [];
  isLoading = false;
  construct: Construct = null;
  history: UserHistory = null;
  exampleFile = 'assets/files/NC_044937.gbk';
  fileList = [];
  endpoint: string;
  view = true;
  search: string;
  text: string = null;

  beforeUpload = (file: File) => {
    this.fileList = [file];
    this.construct = null;
  };

  constructor(
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private router: Router,
    private notify: NotifyService,
    private http: HttpClient
  ) {
    this.endpoint = env.endpoints.api + '/constructs/from-genbank';
  }

  ngOnInit() {
    this.sub = this.route.queryParams.subscribe(params => {
      if (params.specie) this.specie = new Specie();
      this.specie.slug = params.specie || null
    });
    if (this.specie) { this.getSpecie(); }
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

  handleChange({ file }: UploadChangeParam): void {
    if (file.status === 'done') {
      this.construct = new Construct().deserialize(file.response);
    }
  }


  loadExample() {
    this.http.get(this.exampleFile, { observe: 'response', responseType: 'blob' })
      .subscribe((result: HttpResponse<Blob>) => {
        let reader = new FileReader();
        reader.readAsDataURL(result.body);
        reader.onloadend = (ev) => {
          if (ev.target.readyState === window.FileReader.DONE) {
            let data = new Blob([ev.target.result]);
            let arrayOfBlob = new Array<Blob>();
            arrayOfBlob.push(data);
            this.fileList = [new File(arrayOfBlob, result.url.split('/').slice(-1)[0])];
          }
        };
      },
        () => this.notify.warn('Unable to load example file')
      );
  }

  submit() {
    this.response = null;
    this.construct['specie_tax_id'] = this.specie.tax_id;
    this.http
      .post(`${env.endpoints.api}/optimize_seq/from-sketch`, this.construct)
      .subscribe(
        (data: UserHistory) => {
          setTimeout(() => {
            this.router.navigate(['/workspace', data.id]);
          }, 3000);
        },
        err => {
          this.notify.error(err, 'bottom-right');
        }
      );
  }

}
