import { Injectable, Pipe, PipeTransform } from '@angular/core';

import * as moment from 'moment';

const intervals = {
  year: 31536000,
  month: 2592000,
  week: 604800,
  day: 86400,
  hour: 3600,
  minute: 60,
  second: 1
};

@Injectable()
@Pipe({ name: 'date', pure: true })
export class SafeDatePipe implements PipeTransform {
  transform(value: Date | string, ...args: any[]): any {
    if (value) {
      const [format] = args;
      if (format && format.toLowerCase() === 'ago') { return this.dateAgo(value); }
      value = moment(value).format(format);
    }
    return value;
  }

  dateAgo(date: string | number | Date) {
    const seconds = Math.floor((+new Date() - +new Date(date)) / 1000);
    if (seconds < 29) { return 'Just now'; }// less than 30 seconds ago will show as 'Just now'

    let counter: number;
    for (const i of Object.keys(intervals)) {
      counter = Math.floor(seconds / intervals[i]);
      if (counter > 0) {
        return counter === 1 ? `${counter} ${i} ago` : `${counter} ${i}s ago`; // plural or singular
      }
    }
  }
}
