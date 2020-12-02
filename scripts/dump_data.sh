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

 # gunzip -c ~/Data/MetWork/prod_copy/database_backup.dump.gz | psql -v -d metwork -U metwork -h metwork_db

 # SELECT id FROM base_molecule WHERE mol_rdkit="CCCCCCCCCC1=CC(=O)C2=CC=CC=C2N1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1"

 mols = Molecule.objects.raw("DELETE FROM base_molecule WHERE mol_rdkit = 'CCCCCCCCCC1=CC(=O)C2=CC=CC=C2N1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1'")

 pg_dump -v -d metwork -U metwork -t base_molecule > molecule.dump




CCCCCCCCCC1=CC(=O)C2=CC=CC=C2N1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1

CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1

CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1

CCCCCCCCCC1=CC(=O)c2ccccc2[NH]1(C)C(=O)[C@H](Cc1ccccc1)NC(=O)c1ccccc1


CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1

CCCCCCCCCC1=CC(=O)c2ccccc2[NH]1(C)C(=O)[C@H](Cc1ccccc1)NC(=O)c1ccccc1

CCCCCCCCCC1=CC(=O)c2ccccc2[N]1(C)C(=O)[C@H](Cc1ccccc1)NC(=O)c1ccccc1


124414		HPJYSXMDEPPDPV-HKBQPEDESA-N	CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1	CCCCCCCCCC1=CC(=O)c2ccccc2N1(C)C(=O)[C@H](Cc1ccccc1)NC(=O)c1ccccc1
124414		HPJYSXMDEPPDPV-HKBQPEDESA-N	CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1	CCCCCCCCCC1=CC(=O)c2ccccc2N1(C)C(=O)[C@H](Cc1ccccc1)NC(=O)c1ccccc1

CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1

BAD :  CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1
GOOD : CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1
       CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1
       CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[NH]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1


       sed -i -e 's/$BAD_SMILES/$GOOD_SMILES/g' test.txt

       "CCCCCCCCCC1=CC(=O)C2=CC=CC=C2[N+]1(C)C(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CC=CC=C1"