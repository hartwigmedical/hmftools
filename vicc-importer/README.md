# VICC-Importer

VICC is the [Variant Interpretation for Cancer Consortium](https://cancervariants.org). 
One of the aims of this consortium is to create a harmonized meta-knowledgebase from various actual knowledgebases (see [paper](http://dx.doi.org/10.1038/s41588-020-0603-8)).

This module imports the data that has been generated as part of this paper, and in addition does the following:
 *  Determines the event type defined by the [SERVE](../serve/README.md) datamodel.
 *  Provides classification mapping to map the VICC database onto [SERVE](../serve/README.md) actionable events.
 *  Loads up the database into a MySQL database
 
## Loading the VICC database into a MySQL database

In order to load the VICC database into a MySQL database, one needs to:
 * Create the database using create_vicc_database.sql script from the VICC importer resources.
 * Run ViccJsonSQLImporter with the following arguments:
 
Argument  | Description
---|---
vicc_json  | Required: Path to the VICC json file containing the VICC data.
db_url | Required: The URL of the database in which to ingest the data.
db_user | Required: The DB user which has access to the the URL specified
db_pass | Required: The password needed for the DB user to authenticate on the URL.
