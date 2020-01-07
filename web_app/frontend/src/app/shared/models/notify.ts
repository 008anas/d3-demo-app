export class Notify {
    type: NotifyType;
    message: string;
    position: string;
}

export enum NotifyType {
    Success,
    Error,
    Info,
    Warning
}
