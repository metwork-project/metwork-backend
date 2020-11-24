# metwork-backend

This code corespond to the backend part of [MetWork](https://metwork.pharmacie.parisdescartes.fr/) web platform, written in Python with Django framework.

## Conda

Conda is required to run MetWork. If you not ahave already Conda installed, it is recommended to use [miniconda](https://conda.io/miniconda.html).
Python version tested is 3.6.

## Automation scripts

Install conda env and all depedencies :

```
. scripts/install-all.sh 
```

Run development environment :

```
. scripts/run-dev.sh 
```

Run all tests

```
. scripts/run-tests.sh 
```

## Development environement deatil

### Python scripts

The python backend is based on [Django Framework](https://www.djangoproject.com/).The python scripts run via an conda virtual environement (venv). Install the python libraires with  `scripts\install-python-depencies.sh` command.

To run properly, some environment variables have to be set. The easiest way is ti run `. scripts/set-env.sh`. Take a look inside this script to ee how to modifivy environment variables values.

Once venv installed and activated with all dependencies set, run the following :

```
./manage.py [cmd]
```

Refer to Django documentaiton for `[cmd]` values. For example, run `./manage.py runserver` to run the dev server.

### Celery worker

Additional to web server, you must run celery worker to run properly. The easiest way is :

In a new shell if you want to see worker details :

```
. scripts/run-worker.sh 
```

Or in detached mode

```
. scripts/run-worker.sh detached
```

If you want a worker to run tests (pointing to test database) :

```
. scripts/run-worker.sh test
```

### Docker services

Additional services are provided by docker instances :

- Postgres database with RDKIT extension
- Cache
- AMQP Broker

You can run all the services in detached mode with `. scripts/run-docker.sh test` or see inside `docker` folder for further information.
