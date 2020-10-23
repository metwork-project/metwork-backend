FROM ubuntu:16.04

LABEL MAINTAINER="Yann Beauxis <pro@yannbeauxis.net>"

#### ubuntu packages

RUN apt-get update \
  && apt-get install -y \
  # bzip2 locales wget netcat
  wget locales

#### Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

#### User metwork

RUN useradd metwork -m -p 'metwork_pwd'

#### Conda & Python env

RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /root/miniconda.sh \
  && chmod 764 /root/miniconda.sh \
  && /bin/bash /root/miniconda.sh -b -p /opt/conda \
  && rm /root/miniconda.sh \
  && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
  && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
  && rm -R /opt/conda/pkgs/

ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN conda create -n metwork python=3.6

ENV PATH /opt/conda/envs/metwork/bin:$PATH
ENV CONDA_DEFAULT_ENV metwork
ENV CONDA_PREFIX /opt/conda/envs/metwork


COPY conda-requirements.txt /opt/metwork-backend/conda-requirements.txt
COPY pip-requirements.txt /opt/metwork-backend/pip-requirements.txt

RUN echo "conda activate metwork" >> ~/.bashrc \
  && echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/metwork/.bashrc \
  && echo "conda activate metwork" >> /home/metwork/.bashrc

#### Python packages

RUN conda install -c r -c rdkit -c conda-forge -c bioconda --yes --file /opt/metwork-backend/conda-requirements.txt \
  && pip install -r /opt/metwork-backend/pip-requirements.txt \
  && rm -R /opt/conda/pkgs/

#### CFM-ID

RUN wget -q https://metwork.pharmacie.parisdescartes.fr/vendor/cfm_id.tar.gz -O /opt/cfm_id.tar.gz \
  && tar -xf /opt/cfm_id.tar.gz -C /opt/ \
  && rm /opt/cfm_id.tar.gz

RUN mkdir /var/log/metwork/

EXPOSE 80
EXPOSE 5432