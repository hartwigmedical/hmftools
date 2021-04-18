# CKB-Importer

[CKB](https://ckbhome.jax.org) - **C**linical **K**nowledge**b**ase - is a cancer knowledgebase provided by Jackson Lab. 
 
 This module imports the data that is provided in json format (CKB FLEX), and in addition does the following:
  *  For each molecular profile, determines the event type as defined by the [SERVE](../serve/README.md) datamodel.
  *  Provides classification mapping to map all molecular profiles onto [SERVE](../serve/README.md) actionable events.
  *  Loads up the database into a MySQL database
  
 ## Loading the CKB FLEX database into a MySQL database
 
 In order to load the CKB FLEX database into a MySQL database, one needs to:
  * Create the database using generate_ckb_db.sql script from the ckb-importer resources.
  * Run CkbImporterApplication (default class in the ckb-importer jar) with the following arguments:
  
 Argument  | Description
 ---|---
 ckb_dir  | Required: Path to the directory holding the JSON data.
 db_url | Required: The URL of the database in which to ingest the data.
 db_user | Required: The DB user which has access to the the URL specified
 db_pass | Required: The password needed for the DB user to authenticate on the URL.