# metwork-backend

This code corespond to the backend part of [MetWork](https://metwork.pharmacie.parisdescartes.fr/) web platform, written in Python with Django framework.

## Install development environment

### Conda

Conda is required to run MetWork. If you not ahave already Conda installed, it is recommended to use [miniconda](https://conda.io/miniconda.html).
Python version tested is 3.6.

You must create a dedicated conda env for metwork :

```
conda create -n metwork python=3.6
conda active metwork
```

### Dependencies

From repository root folder :

```
conda activate metwork
bash install/python-dependencies.sh 
bash install/cfm_id.sh 
```

### Docker images

You must install [docker](https://docs.docker.com/engine/install/) and [docker-compose](https://docs.docker.com/compose/install/). 

Then from `docker` folder :

```
docker-compose pull
```

### Init db

from root folder :

```
cd docker
docker-compose up database -d
cd ..
bash set-env.sh && ./manage.py migrate
cd docker
docker-compose stop database
```


## Run dev environment

1. Docker instances

Then from `docker` folder :

```
docker-compose up 
```

2. Worker

From repository root folder.

```
conda activate metwork
bash run-worker.sh
```

3. API

From repository root folder.

```
conda activate metwork
bash run-api.sh
```
