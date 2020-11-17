#!/bin/bash
set -e

psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" <<-EOSQL
  ALTER ROLE metwork SET client_encoding TO 'utf8';
  ALTER ROLE metwork SET default_transaction_isolation TO 'read committed';
  ALTER ROLE metwork SET timezone TO 'UTC';
  ALTER ROLE metwork SUPERUSER;
  GRANT ALL PRIVILEGES ON DATABASE metwork TO metwork;
EOSQL

