// This file can be replaced during build by using the `fileReplacements` array.
// `ng build --prod` replaces `environment.ts` with `environment.prod.ts`.
// The list of file replacements can be found in `angular.json`.

export const environment = {
  production: false,
  email: 'bioroboost@crg.es',
  name: 'SQRUTINY',
  colors: {
    main: '#007bff'
  },
  sentry: {
    dsn: 'https://f35c3be6cd794b069c754d34e0daa7a9@sentry.io/1552997'
  },
  routes: {
    optimize: {
      root: 'optimize',
      sketcher: 'sketcher',
      from_file: 'file'
    },
    construct: {
      root: 'construct'
    },
    workspace: {
      root: 'workspace',
      detail: ':id'
    },
    documentation: 'documentation',
    home: 'home',
    contact: 'contact',
    error404: '404',
    error500: '500'
  },
  endpoints: {
    api: 'http://localhost:8000/api/v1',
  }
};

/*
 * For easier debugging in development mode, you can import the following file
 * to ignore zone related error stack frames such as `zone.run`, `zoneDelegate.invokeTask`.
 *
 * This import should be commented out in production mode because it will have a negative impact
 * on performance if an error is thrown.
 */
import 'zone.js/dist/zone-error';  // Included with Angular CLI.
