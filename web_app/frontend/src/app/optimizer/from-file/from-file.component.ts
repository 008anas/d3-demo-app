import { Component, OnInit, OnDestroy } from '@angular/core';
import { Subscription } from 'rxjs/internal/Subscription';
import { ActivatedRoute, Router } from '@angular/router';
import { HttpClient, HttpResponse, HttpEvent, HttpEventType, HttpRequest } from '@angular/common/http';
import { finalize } from 'rxjs/operators';

import { UploadXHRArgs } from 'ng-zorro-antd/upload';
import { NzModalService } from 'ng-zorro-antd/modal';
import { NzMessageService } from 'ng-zorro-antd/message';

import { environment as env } from '@env/environment';
import { Specie } from '@models/specie';
import { Construct } from '@models/construct';
import { SpecieService } from '@services/specie.service';
import { UserHistory } from 'app/workspace/shared/user-history';
import { TextModalComponent } from '@components/text-modal/text-modal.component';
import { SqrutinyService } from '@services/sqrutiny.service';
import { LoaderService } from '@services/loader.service';
import { NavService } from '@services/nav.service';

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
  isSubmitting = false;

  customReq = (item: UploadXHRArgs) => {
    this.loader.startLoading();
    this.fileList = [item.file];
    const formData = new FormData();
    formData.append('file', item.file as any);
    const req = new HttpRequest('POST', item.action!, formData, {
      reportProgress: true,
      withCredentials: true
    });
    this.response = null;
    this.construct = null;
    return this.http.request(req)
      .pipe(finalize(() => this.loader.stopLoading()))
      .subscribe(
        (event: HttpEvent<any>) => {
          if (event.type === HttpEventType.UploadProgress) {
            if (event.total! > 0) {
              (event as any).percent = (event.loaded / event.total!) * 100;
            }
            item.onProgress!(event, item.file!);
          } else if (event instanceof HttpResponse) {
            item.file.status = 'done';
            item.onSuccess!(event.body, item.file!, event);
            this.construct = new Construct().deserialize(event.body);
            this.notify.success('GenBank file parsed successfully!');
          }
        },
        err => this.response = err
      );
  }

  constructor(
    private route: ActivatedRoute,
    private specieSrvc: SpecieService,
    private sqrutinySrvc: SqrutinyService,
    private router: Router,
    private notify: NzMessageService,
    private http: HttpClient,
    private modal: NzModalService,
    private loader: LoaderService,
    private navSrvc: NavService
  ) {
    this.endpoint = env.endpoints.api + '/constructs/from-genbank';
  }

  ngOnInit() {
    this.sub = this.route.queryParams.subscribe(params => {
      if (params.specie) { this.specie = new Specie(); }
      this.specie.slug = params.specie || null;
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


  loadExample() {
    this.http.get(this.exampleFile, { observe: 'response', responseType: 'blob' })
      .subscribe((result: HttpResponse<Blob>) => {
        const reader = new FileReader();
        reader.readAsDataURL(result.body);
        reader.onloadend = (ev) => {
          if (ev.target.readyState === window.FileReader.DONE) {
            const data = new Blob([ev.target.result]);
            const arrayOfBlob = new Array<Blob>();
            arrayOfBlob.push(data);
            this.fileList = [new File(arrayOfBlob, result.url.split('/').slice(-1)[0])];
          }
        };
      },
        () => this.notify.warning('Unable to load example file')
      );
  }

  submit() {
    this.response = null;
    this.isSubmitting = true;
    this.construct.specie_tax_id = this.specie.tax_id;
    this.sqrutinySrvc.fromSketch(this.construct)
      .subscribe(
        (data: UserHistory) => {
          this.navSrvc.updateBadge();
          this.notify.loading('Your job has been submitted! In a moment you will be redirected. Please be patient');
          setTimeout(() => {
            this.router.navigate(['/workspace', data.id]);
            this.isSubmitting = false;
          }, 3000);
        },
        err => {
          this.isSubmitting = false;
          this.notify.error(err);
        }
      );
  }

  textModal(str: string) {
    this.modal.create({
      nzContent: TextModalComponent,
      nzWrapClassName: 'center-modal',
      nzComponentParams: {
        txt: str
      },
      nzFooter: null
    });
  }

}
