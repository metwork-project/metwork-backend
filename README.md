# metwork-backend

This code corespond to the backend part of [MetWork](https://metwork.pharmacie.parisdescartes.fr/) web platform, written in Python with Django framework.

## Python environment

It's recommended to use python dedicated environment with tool such as [miniconda](https://conda.io/miniconda.html). Python version tested is 3.6.4.

## Dependencies

Python libraries used are specify in conda-requirements.txt for those to be installed by conda and pip-requirements.txt for those to be installed by pip. The [django_rdkit](https://django-rdkit.readthedocs.io/en/latest/) app is also used.

[ChemAxon Reactor](https://chemaxon.com/products/reactor) is used and has to be installed. You need to procure your own license and copy it to the path of the  CHEMAXON_LICENSE_URL environment variable.

[CFM-ID](http://cfmid.wishartlab.com/) is used and a [compiled version](https://metwork.pharmacie.parisdescartes.fr/vendor/cfm_id.tar.gz) as to be extract and copy to the path of the METWORK_CFM_ID_PATH environment variable.

## Environment variables

The following environment variables have to be set :

- METWORK_BACKEND_PATH=/opt/metwork-backend
- METWORK_DB_PASSWORD=METWORK_DB_PASSWORD
- METWORK_BROKER_PASSWORD=METWORK_BROKER_PASSWORD
- METWORK_SECRET_KEY=METWORK_SECRET_KEY
- METWORK_ALLOWED_HOSTS=0.0.0.0
- METWORK_DATA_FILES_PATH=/srv/metwork/files
- METWORK_CFM_ID_PATH=/opt/cfm_id/
- CHEMAXON_LICENSE_URL=/opt/chemaxon/license/license.cxl

## Services required

This app need to run :

- A postgreSQL database with the [RDKit cartridge](http://www.rdkit.org/docs/Cartridge.html)
- A Memcached service
- A Rabbitmq service
- Celery workers connected to the queues of the CELERY_QUEUES django settings
