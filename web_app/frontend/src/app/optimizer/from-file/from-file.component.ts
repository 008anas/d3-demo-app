import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';
import { ActivatedRoute, Router } from '@angular/router';
import { HttpClient, HttpEvent, HttpEventType, HttpRequest, HttpResponse } from '@angular/common/http';
import { finalize } from 'rxjs/operators';

import { UploadChangeParam } from 'ng-zorro-antd/upload';
import { UploadXHRArgs } from 'ng-zorro-antd/upload';

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
  specie: Specie = null;
  species: Specie[] = [];
  isLoading = false;
  history: UserHistory = null;
  isUploading = false;
  construct: Construct = null;
  exampleFile = 'assets/files/NC_044937.gbk';
  endpoint: string;
  view = true;
  search: string;
  isProcessing = false;

  customReq = (item: UploadXHRArgs) => {
    const formData = new FormData();
    // tslint:disable-next-line:no-any
    formData.append('file', item.file as any);

    const req = new HttpRequest('POST', item.action!, formData, {
      reportProgress: true,
      withCredentials: true
    });
    // Always returns a `Subscription` object. nz-upload would automatically unsubscribe it at correct time.
    return this.http.request(req).subscribe(
      (event: HttpEvent<any>) => {
        switch (event.type) {
          case HttpEventType.UploadProgress:
            if (event.total! > 0) {
              // tslint:disable-next-line:no-any
              (event as any).percent = (event.loaded / event.total!) * 100;
            }
            item.onProgress!(event, item.file!);
            break;
          case HttpEventType.Response:
            this.construct = new Construct().deserialize(event.body);
            item.onSuccess!(event.body, item.file!, event);
            break;
        }
      },
      err => console.log(err)
    );
  };

  constructor(
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private router: Router,
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
    const status = file.status;
    if (status !== 'uploading') {
      this.isProcessing = true;
    }
    if (status === 'done') {
      this.isProcessing = false;
    }
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
