BACKUP_DATE=$(date +"%Y-%m-%d")

# pg_dump -v -h metwork_db -d metwork -U metwork > database_backup_$BACKUP_DATE.dump

DATA_FOLDER=/srv/metwork
DB_FOLDER=$DATA_FOLDER/database_backup
DB_FILE=$DB_FOLDER/database_backup_$BACKUP_DATE.dump.gz

pg_dump -v -d metwork -U metwork | gzip > $DB_FILE

ln -s $DB_FILE $DB_FOLDER/database_backup.dump.gz

FILES_FOLDER=$DATA_FOLDER/files_backup
FILES_FILE=$DB_FOLDER/files_backup_$BACKUP_DATE.tar.gz

tar -zcvf $FILES_FILE $DATA_FOLDER/files

ln -s $FILES_FILE $FILES_FOLDER/database_backup.dump