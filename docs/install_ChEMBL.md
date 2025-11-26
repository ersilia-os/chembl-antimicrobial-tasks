# Instructions to install the ChEMBL database in a local computer

These instructions are valid for Ubuntu Linux.

1. Install postgreSQL database system.
```
$ sudo apt update
$ sudo apt install postgresql postgresql-contrib
$ sudo apt install postgresql-client
```

2. Start service.

(Note: **you will have to run this again if you restart your computer**, or alternatively have it run automatically on startup)
```
$ sudo service postgresql start
```

The user "postgres" is created automatically and is the DB administrator. To do any configuration, you need to run as user postgres (sudo -u postgres).

3. Create database for ChEMBL (make sure to include the final semicolon in all SQL commands)
```
$ sudo -u postgres psql
postgres=# create database chembl_36;
postgres=# \q
```
At that step, you might encounter issues due to the permissions in your user. Please look at this [thread](https://dba.stackexchange.com/questions/54242/cannot-start-psql-with-user-postgres-could-not-change-directory-to-home-user) or this [thread](https://www.reddit.com/r/PostgreSQL/comments/vc2nbw/ubuntukubuntu_2204_permission_denied_how_to_solve/) to learn more and solve them.

4. Download the ChEMBL data in postgres format
From this site: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
Download the file: `chembl_36_postgresql.tar.gz`
Decompress it to get the file `chembl_36_postgresql.dmp`

When decompressing:
````
tar -xvzf chembl_36_postgresql.tar.gz
````
the command will create a new folder in your working directory `chembl_36/`.

This directory contains the database dump and installation notes provided by ChEMBL. You only need this folder temporarily to run `pg_restore` (see next step). After importing the database into PostgreSQL, the `chembl_36/` directory may be deleted safely.

5. Load database contents using the downloaded file. NOTE: The database requires about 23 GB of disk space.

This step will take time: about 15-25 minutes on an SSD disk, considerably longer on a traditional HDD disk. It took roughly 20 min in Ersilia's machine `herbert`.  
```
$ sudo -u postgres pg_restore --no-owner --verbose -U postgres -d chembl_36 chembl_36_postgresql.dmp
```

6. Test: Enter the database for querying. Run a test query. You should get the list of assay types.
```
$ sudo -u postgres psql chembl_36
chembl_36=# select * from assay_type;
 assay_type |   assay_desc    
------------+-----------------
 A          | ADME
 B          | Binding
 F          | Functional
 P          | Physicochemical
 T          | Toxicity
 U          | Unassigned
(6 rows)

chembl_36=# \q
```

So far the database is installed and it works with user postgres. For security, we prefer to create a specific DB user
with limited rights, called "chembl_user", to use in our programs. The password for this DB user
will be "aaa".

7. Create a specific user for querying this database
```
$ sudo -u postgres createuser --interactive -P
Enter name of role to add: chembl_user
Enter password for new role: aaa
Enter it again: aaa
Shall the new role be a superuser? (y/n) n
Shall the new role be allowed to create databases? (y/n) n
Shall the new role be allowed to create more new roles? (y/n) n
```

8. Grant access to the chembl_36 tables to user chembl_user
```
$ sudo -u postgres psql chembl_36
chembl_36=# GRANT SELECT ON ALL TABLES IN SCHEMA public TO chembl_user;
chembl_36=> \q
```

9. Test: Connect to database as user chembl_user. Check the structure of the database, which is explained in detail in [ChEMBL's site](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.txt).
```
$ psql -h localhost -p 5432 -U chembl_user chembl_36
Password for user chembl_user: (use the password you created before)
psql (12.13 (Ubuntu 12.13-0ubuntu0.20.04.1))
SSL connection (protocol: TLSv1.3, cipher: TLS_AES_256_GCM_SHA384, bits: 256, compression: off)
Type "help" for help.
chembl_36=> \d
                    List of relations
 Schema |             Name             | Type  |  Owner   
--------+------------------------------+-------+----------
 public | action_type                  | table | postgres
 public | activities                   | table | postgres
 public | activity_properties          | table | postgres
 public | activity_smid                | table | postgres
 public | activity_stds_lookup         | table | postgres
 public | activity_supp                | table | postgres
 public | activity_supp_map            | table | postgres
 public | assay_class_map              | table | postgres
 public | assay_classification         | table | postgres
 public | assay_parameters             | table | postgres
 public | assay_type                   | table | postgres
 public | assays                       | table | postgres
 public | atc_classification           | table | postgres
 public | binding_sites                | table | postgres
 public | bio_component_sequences      | table | postgres
 public | bioassay_ontology            | table | postgres

 ...

chembl_36=> \q
```

**DONE!**

