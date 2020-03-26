"""
Django settings for sqrutiny project.

Generated by 'django-admin startproject' using Django 2.2.6.
"""

import os
from configparser import RawConfigParser

import sentry_sdk
from sentry_sdk.integrations.django import DjangoIntegration

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

config = RawConfigParser()

config.read(os.path.join(BASE_DIR, '.env.ini'))

# SECURITY WARNING: keep the secret key used in production secret!
if config.has_option('MAIN', 'SECRET_KEY'):
    SECRET_KEY = config.get('MAIN', 'SECRET_KEY')
else:
    SECRET_KEY = 'mncl3yc&l&a4@6c_)9*n!tn+p@8o7zs^i19w38zfsk6t%f0r65'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = (config.get('MAIN', 'DEBUG').lower() in ['true', '1', 'yes'])

APP_URL = config.get('MAIN', 'APP_URL')

APPEND_SLASH = False

# EMAIL
EMAIL_HOST_USER = 'sqrutiny@crg.es'
EMAIL_HOST_PASSWORD = ''

LOGGING_DIR = os.path.join(BASE_DIR, 'logs')

ALLOWED_HOSTS = ['*']

# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'corsheaders',
    'django_rq',
    'rest_framework',
    'core',
    'rest_framework_tracking',
    'app.specie.apps.SpeciesConfig',
    'app.genetic_element.apps.GeneticElementsConfig',
    'app.contact.apps.ContactConfig',
    'app.construct.apps.ConstructConfig',
    'app.workspace.apps.WorkspaceConfig',
    'app.parameter.apps.ParameterConfig'
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'corsheaders.middleware.CorsMiddleware'
]

# REST FRAMEWORK configuration
REST_FRAMEWORK = {
    'DEFAULT_RENDERER_CLASSES': (
        'rest_framework.renderers.JSONRenderer',
        'rest_framework.renderers.BrowsableAPIRenderer'
    ),
    'DEFAULT_THROTTLE_CLASSES': (
        'rest_framework.throttling.AnonRateThrottle',
        'rest_framework.throttling.UserRateThrottle'
    ),
    'DEFAULT_THROTTLE_RATES': {
        'anon': '10000/hour',
        'user': '10000/hour',
    }
}

# CORS configuration
CORS_ORIGIN_ALLOW_ALL = True
CORS_ALLOW_CREDENTIALS = True

SESSION_ENGINE = 'django.contrib.sessions.backends.cached_db'
SESSION_EXPIRE_AT_BROWSER_CLOSE = False

ROOT_URLCONF = 'sqrutiny.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'sqrutiny.wsgi.application'

# Database

if config.has_section('DATABASE'):
    DATABASES = {
        'default': {
            'ENGINE': config.get('DATABASE', 'DB_ENGINE'),
            'NAME': config.get('DATABASE', 'DB_NAME'),
            'USER': config.get('DATABASE', 'DB_USER'),
            'PASSWORD': config.get('DATABASE', 'DB_PASS'),
            'HOST': config.get('DATABASE', 'DB_HOST'),
            'PORT': config.get('DATABASE', 'DB_PORT')
        }
    }
else:
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': os.path.join(BASE_DIR, 'db.sqlite3')
        }
    }

# Password validation

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# Static files (CSS, JavaScript, Images)

STATIC_URL = '/static/'

MEDIA_ROOT = os.path.join(BASE_DIR, 'media')
MEDIA_URL = '/media/'

# RQ WORKER

RQ_QUEUES = {
    'default': {
        'URL': 'redis://localhost:6379',
        # 'PASSWORD': 'gJKUadGMAaSUeYEKuDSEX0nHPjzfC3DY1kakCnmX3XmkVNboqIVqT1+jw35SHGcvb0efIqYwb9Apm32h',
        'PORT': 6379,
        'DB': 0,
        'DEFAULT_TIMEOUT': 360,
    }
}

RQ_SHOW_ADMIN_LINK = True

sentry_sdk.init(
    dsn='https://18fa303e3ac248acac6faa0455ab6160@sentry.io/5174972',
    integrations=[DjangoIntegration()],
    send_default_pii=True
)
