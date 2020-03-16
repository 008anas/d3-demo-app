import { Injectable } from '@angular/core';

import { saveAs } from 'file-saver';
import * as XLSX from 'xlsx';

const EXCEL_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
const EXCEL_EXTENSION = '.xlsx';

@Injectable({
  providedIn: 'root'
})
export class FileService {

  constructor() { }

  public exportAsExcelFile(json: any[], excelFileName: string): void {
    const worksheet: XLSX.WorkSheet = XLSX.utils.json_to_sheet(json);
    const workbook: XLSX.WorkBook = { Sheets: { data: worksheet }, SheetNames: ['data'] };
    const excelBuffer: any = XLSX.write(workbook, { bookType: 'xlsx', type: 'array' }); // type: 'buffer'
    this.saveFileAs(excelBuffer, EXCEL_TYPE, excelFileName + '_export_' + new Date().getTime() + EXCEL_EXTENSION);
  }

  public saveFileAs(buffer: any, type: string, fileName: string): void {
    const data: Blob = new Blob([buffer], {type});
    saveAs(data, fileName);
  }
}
