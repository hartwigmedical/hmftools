# HMF Patients Database

## Create Database

The following commands will create a user with write permissions, a user with read permissions and a database.  

```
mysql> ​CREATE USER 'purple_writer'@'localhost' IDENTIFIED WITH mysql_native_password BY 'purple_writer_password'; 
Query OK, 0 rows affected (0.00 sec)
mysql> ​CREATE USER 'purple_reader'@'localhost' IDENTIFIED WITH mysql_native_password BY 'purple_reader_password'; 
Query OK, 0 rows affected (0.00 sec)
mysql> CREATE DATABASE patientdb; 
Query OK, 1 row affected (0.00 sec)
mysql> GRANT ALL on patientdb.* TO 'purple_writer'@'localhost'; 
Query OK, 0 rows affected (0.00 sec)
mysql> GRANT SELECT on patientdb.* TO 'purple_reader'@'localhost'; 
Query OK, 0 rows affected (0.00 sec)
```

## Create Tables
If creating a database from scratch, execute the [generate_database.sql](../patient-db/src/main/resources/generate_database.sql) script from the command with the following. 
Note that you will be prompted for a password:

```
mysql -u writer -p < generate_database.sql
```


### Update Tables
Updates to the schema are recorded in the [patches directory](../patient-db/src/main/resources/patches/patientdb). 
They can be executed in a similar manner as above.


## Data loaders

Data will be deleted before new records are inserted. The loaders do not support updating records.

### Purple data loading

```
java -cp patient-db.jar com.hartwig.hmftools.patientdb.LoadPurpleData \ 
    -sample COLO829T \
    -reference COLO829R \
    -purple_dir /path/COLO829/purple \
    -db_user writer -db_pass writer_password \
    -db_url mysql://localhost:3306/patientdb?serverTimezone=UTC
```

Either somatic or germline data only can be loaded if the configs -somatic_only and -germline_only are included.

Note that if the somatic variants also contain a reference or an rna sample (or both) these can be loaded by supplying the optional arguments `reference` and `rna` respectively, eg:


### Linx data loading

```
java -cp patient-db.jar com.hartwig.hmftools.patientdb.LoadLinxData \ 
    -sample COLO829T \
    -linx_dir /path/COLO829/linx \
    -data_type [both=default, somatic, germline] \
    -db_user writer -db_pass writer_password \
    -db_url mysql://localhost:3306/patientdb?serverTimezone=UTC
```

Either somatic or germline data only can be loaded if the data_type config is provided.

## Use Queries

After connecting to the PURPLE database (eg `mysql -u purple_reader -p -d patientdb`), you can query the following PURPLE tables:

```
SELECT * FROM purity WHERE sampleId = 'COLO829T';
SELECT * FROM purityRange WHERE sampleId = 'COLO829T';
SELECT * FROM copyNumber WHERE sampleId = 'COLO829T';
SELECT * FROM copyNumberGermline WHERE sampleId = 'COLO829T';
SELECT * FROM geneCopyNumber WHERE sampleId = 'COLO829T';
SELECT * FROM structuralVariant WHERE sampleId = 'COLO829T';
SELECT * FROM somaticVariant WHERE sampleId = 'COLO829T';
```

With all the resources of SQL available, it is possible to construct powerful and informative queries such as this example of how one might 
calculate the microsatellite status of a sample by examining its indels:

```
SELECT sampleId, count(*)/2859 as indelsPerMb, if(count(*)/2859 > 4, "MSI", "MSS" ) AS status 
FROM somaticVariant
WHERE filter = 'PASS'
 AND type = 'INDEL' AND repeatCount >= 4 AND length(alt) <= 50 AND length(ref) <= 50
 AND ((length(repeatSequence) BETWEEN 2 AND 4 ) OR
	 (length(repeatSequence) = 1 AND repeatCount >= 5))
 AND sampleId IN ('COLO829T')
GROUP BY sampleId
ORDER BY 2 DESC;
```  

Similarly, it is possible to query the tumor mutation burden of all samples with a query such as:

```
SELECT sampleId, count(*)/2859 AS TMB, IF (count(*)/2859 > 10, "High", "Low") AS status
FROM somaticVariant 
WHERE filter = 'PASS'
GROUP BY 1;
```

Regions of kataegis within a sample are queried with:
```
SELECT kataegis, min(chromosome) as chromosome, min(position) as start, max(position) as end,  
       count(*), round((max(position) - min(position)) / (count(*) - 1))  as avgDistance
FROM somaticVariant 
WHERE sampleId = 'COLO829T' AND kataegis <> ''
GROUP BY kataegis
ORDER BY count(*) DESC;
```

