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
 
 Resources for v1.0 can be found on https://github.com/hartwigmedical/hmftools/releases/tag/ckb-importer-v1.0

## Version History and Download Links
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/ckb-importer-v1.3)
  - Implemention treatment approches of the evidence into the datamodel
  - Study phase could be nullable 
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/ckb-importer-v1.2)
  - Support for java11
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/ckb-importer-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/ckb-importer-v1.0)
- Initial release. 