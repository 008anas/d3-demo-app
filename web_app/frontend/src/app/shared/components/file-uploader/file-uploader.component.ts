import { Component, Input, Output, EventEmitter } from '@angular/core';
import { HttpClient, HttpEventType } from '@angular/common/http';

@Component({
  selector: 'sqy-file-uploader',
  templateUrl: './file-uploader.component.html',
  styleUrls: ['./file-uploader.component.scss']
})
export class FileUploaderComponent {

  @Input() multiple = false;
  @Input() endpoint: string;
  @Output() onUpload: EventEmitter<any> = new EventEmitter<any>();
  response: any = null;
  files: any[] = [];

  constructor(private http: HttpClient) { }

  /**
   * on file drop handler
   */
  onFileDropped($event) {
    this.prepareFilesList($event);
  }

  /**
   * handle file from browsing
   */
  fileBrowseHandler(files) {
    this.prepareFilesList(files);
  }

  /**
   * Delete file from files list
   * @param index (File index)
   */
  deleteFile(index: number) {
    this.files.splice(index, 1);
  }

  /**
   * Simulate the upload process
   */
  uploadFilesSimulator(index: number) {
    if (!this.files[index]) { return; }
    const params = new FormData();
    params.append('file', this.files[index]);
    this.response = null;
    this.http.post(this.endpoint, params, { reportProgress: true, observe: 'events' })
      .subscribe(
        event => {
          switch (event.type) {
            case HttpEventType.Sent:
              console.log('Request has been made!');
              break;
            case HttpEventType.ResponseHeader:
              console.log('Response header has been received!');
              break;
            case HttpEventType.UploadProgress:
              this.files[index].progress = Math.round(event.loaded / event.total * 100);
              break;
            case HttpEventType.Response:
              this.onUpload.emit(event.body);
              break;
          }
        },
        err => this.response = err
      );
  }

  /**
   * Convert Files list to normal array list
   * @param files (Files List)
   */
  prepareFilesList(files: Array<any>) {
    for (const item of files) {
      item.progress = 0;
      this.files.push(item);
      this.uploadFilesSimulator(this.files.length - 1);
    }
  }

}
