export const environment = {
  production: true,
  email: 'bioroboost@crg.es',
  colors: {
    main: '#007bff'
  },
  sentry: {
    dsn: 'https://f35c3be6cd794b069c754d34e0daa7a9@sentry.io/1552997'
  },
  routes: {
    optimize: {
      root: 'optimize/:specie',
      sketcher: 'sketcher',
      sketcher_with_id: 'sketcher/:uuid'
    },
    construct: {
      root: 'construct'
    },
    workspace: {
      root: 'workspace',
      detail: ':uuid'
    },
    documentation: 'documentation',
    home: 'home',
    contact: 'contact',
    error404: '404',
    error500: '500'
  },
  endpoints: {
    api: 'http://bioroboost.crg.es/api/v1',
  }
};
