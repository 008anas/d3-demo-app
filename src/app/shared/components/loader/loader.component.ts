import { Component, AfterViewInit, OnDestroy, ElementRef, ChangeDetectorRef } from '@angular/core';
import { LoaderService } from '@services/loader.service';
import { Subscription } from 'rxjs/internal/Subscription';
import { debounceTime } from 'rxjs/operators';

@Component({
  selector: 'sqy-loader',
  templateUrl: './loader.component.html',
  styleUrls: ['./loader.component.scss']
})
export class LoaderComponent implements AfterViewInit, OnDestroy {

  debounceTime = 200;
  loading = false;
  loadingSubscription: Subscription;

  constructor(
    private loaderSrvc: LoaderService,
    private elmRef: ElementRef,
    private changeDetectorRef: ChangeDetectorRef
  ) { }

  ngAfterViewInit(): void {
    this.elmRef.nativeElement.style.display = 'none';
    this.loadingSubscription = this.loaderSrvc.status()
      .pipe(debounceTime(this.debounceTime))
      .subscribe(
        (status: boolean) => {
          this.elmRef.nativeElement.style.display = status ? 'block' : 'none';
          this.changeDetectorRef.detectChanges();
        }
      );
  }

  ngOnDestroy() {
    this.loadingSubscription.unsubscribe();
  }

}
