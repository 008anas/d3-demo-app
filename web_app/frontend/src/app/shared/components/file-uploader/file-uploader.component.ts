import { Component, OnInit, Input } from '@angular/core';

@Component({
  selector: 'sqy-file-uploader',
  templateUrl: './file-uploader.component.html',
  styleUrls: ['./file-uploader.component.scss']
})
export class FileUploaderComponent implements OnInit {

  @Input() multiple =  false;

  constructor() { }

  ngOnInit() {
  }

}
