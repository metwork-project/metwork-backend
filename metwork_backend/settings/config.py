# coding=utf-8
import os
import configparser
from .utils import get_env

config = configparser.RawConfigParser()
config.read(get_env('METWORK_CONFIG_PATH'))

METWORK_CONF = config