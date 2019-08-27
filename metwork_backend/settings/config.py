# coding=utf-8
import os
import configparser

config = configparser.RawConfigParser()
config.read(os.environ['METWORK_CONFIG_PATH'])

METWORK_CONF = config