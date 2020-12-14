# Install worker

{WORKER_NUMER} has to be replaced with the integer number of worker, i.e. `2` for the worker with the local IP `192.168.0.12`

## Ubuntu server

Current version of ubuntu server used is 20.04.

Hostname : metwork-worker-{WORKER_NUMBER}
User : metwork

Use LVM without encryption.

On install configure network only for eno2 and local network (192.168.0.0).

## Network configuration

Edit configuration file : 

```
sudo nano /etc/netplan/00-installer-config.yaml
```

```yaml
--8<-- "install/00-installer-config.yaml"
```

## Update packages

```
sudo apt-get update && sudo apt-get upgrade
sudo apt-get install -y sysstat nfs-common supervisor
```

## SSH

From worker :

```
ssh-keygen
# default location : /home/metwork/.ssh/id_rsa
# No passphrase
```

From metwork-main :

```
ssh-copy-id 192.168.0.1{WORKER_NUMBER}

```

## Monitoring

### Configure Systat


```
sudo nano /etc/default/sysstat
```

``` hl_lines="9"
#
# Default settings for /etc/init.d/sysstat, /etc/cron.d/sysstat
# and /etc/cron.daily/sysstat files
#

# Should sadc collect system activity informations? Valid values
# are "true" and "false". Please do not put other values, they
# will be overwritten by debconf!
ENABLED="true"
```

### Modify daily cron for systat

```
sudo nano /etc/cron.daily/sysstat
```

```bash hl_lines="17 19 21 23"
--8<-- "install/cron.sysstat"
```

Then (if mandatory ??) :

```
sudo service cron reload
```

## NFS

Add the following line to `/etc/fstab` :

```
192.168.0.1:/srv/metwork/files/ /srv/metwork/files/ nfs defaults,user,auto,noatime,intr 0 0 
192.168.0.1:/opt/metwork-backend/ /opt/metwork-backend/ nfs defaults,user,auto,noatime,intr 0 0
192.168.0.1:/srv/metwork/conf/ /srv/metwork/conf/ nfs defaults,user,auto,noatime,intr 0 0
192.168.0.1:/opt/cfm_id/ /opt/cfm_id/ nfs defaults,user,auto,noatime,intr 0 0
```

Then :

```
sudo mkdir /srv/metwork
sudo mkdir /srv/metwork/files
sudo mkdir /srv/metwork/conf
sudo chown -R metwork:metwork /srv/metwork
sudo mkdir /opt/cfm_id
sudo chown -R metwork:metwork /opt/cfm_id/
sudo mkdir /opt/metwork-backend && sudo chown metwork:metwork /opt/metwork-backend
sudo mount -a
```

## Python scripts

```
CONDA_FILE=Miniconda3-py38_4.9.2-Linux-x86_64.sh
wget https://repo.anaconda.com/miniconda/$CONDA_FILE $HOME
bash $HOME/$CONDA_FILE # Use default config
rm $CONDA_FILE
```

Log out and in to activate conda

```
conda config --set auto_activate_base false
conda deactivate
conda create -y -n metwork python=3.6
conda activate metwork
sudo mkdir /etc/metwork
sudo wget https://raw.githubusercontent.com/metwork-project/metwork-backend/scale-up-preparation/conda-requirements.txt -O /etc/metwork/conda-requirements.txt
sudo wget https://raw.githubusercontent.com/metwork-project/metwork-backend/scale-up-preparation/pip-requirements.txt -O /etc/metwork/pip-requirements.txt
conda install -c r -c rdkit -c conda-forge -c bioconda --yes --file /etc/metwork/conda-requirements.txt
pip install -r /etc/metwork/pip-requirements.txt
```

```
sudo nano /etc/metwork/run_worker.sh
```

```
#! bash
METWORK_VERSION=$(cat /opt/metwork-backend/VERSION)
/home/metwork/miniconda3/envs/metwork/bin/celery worker -A metwork_backend -Q $1.$METWORK_VERSION -l INFO
```

```
sudo nano /etc/supervisor/conf.d/metwork_worker.conf
```

    --8<-- "install/metwork_worker.run.conf"

```
sudo mkdir /var/log/metwork
sudo service supervisor restart
sudo supervisorctl status metwork_worker # Check if running
```

Check on [http://194.254.93.160:15672/#/queues/metwork/run.0.4.4](http://194.254.93.160:15672/#/queues/metwork/run.0.4.4) if worker running.

> Don't forget to release internet connection when finished
