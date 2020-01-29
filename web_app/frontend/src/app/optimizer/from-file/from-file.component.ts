import { Component, OnInit } from '@angular/core';

import { HttpClient, HttpEventType } from '@angular/common/http';
import { FormGroup, FormBuilder } from '@angular/forms';

@Component({
  selector: 'sqy-from-file',
  templateUrl: './from-file.component.html',
  styleUrls: ['./from-file.component.scss']
})
export class FromFileComponent implements OnInit {

  progress = 0;
  form: FormGroup;
  response: any = null;

  constructor(private http: HttpClient, private formBuilder: FormBuilder) { }

  ngOnInit() {
    this.form = this.formBuilder.group({
      file: ['']
    });
  }

  onChange(event: { target: { files: any[]; }; }) {
    if (event.target.files.length > 0) {
      const file = event.target.files[0];
      this.form.get('file').setValue(file);
    }
  }

  onSubmit() {
    const params = new FormData();
    params.append('file', this.form.get('file').value);
    this.response = null;
    this.progress = 0;
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
              this.progress = Math.round(event.loaded / event.total * 100);
              console.log(`Uploaded! ${this.progress}%`);
              break;
            case HttpEventType.Response:
              console.log('User successfully created!', event.body);
          }
        },
        err => {
          this.response = err;
          console.log(err)
        }
      )
  }

}
