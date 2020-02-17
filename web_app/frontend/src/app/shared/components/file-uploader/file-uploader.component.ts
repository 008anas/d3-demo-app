import { Component, Input } from '@angular/core';
import { HttpClient, HttpEventType } from '@angular/common/http';

@Component({
  selector: 'sqy-file-uploader',
  templateUrl: './file-uploader.component.html',
  styleUrls: ['./file-uploader.component.scss']
})
export class FileUploaderComponent {

  @Input() multiple = false;

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
    if (index === this.files.length) {
      return;
    } else {
      const params = new FormData();
      params.append('file', this.files[index].value);
      this.http.post('http://127.0.0.1:8000/api/v1/optimize_seq/from-file', params, { reportProgress: true, observe: 'events' })
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
                break;
            }
          },
          err => {
            // this.response = err.msg;
            console.log(err);
          }
        );
    }
  }

  /**
   * Convert Files list to normal array list
   * @param files (Files List)
   */
  prepareFilesList(files: Array<any>) {
    for (const item of files) {
      item.progress = 0;
      this.files.push(item);
    }
    this.uploadFilesSimulator(0);
  }

}
