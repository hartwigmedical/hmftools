# SV Analyser

## Introduction

This tool will load Structural Variants from the DB and perform a collection of analytical routines on the data.

The output of this analysis will then be written to file, and loaded into new DB tables.

## Resources

None

## Dependencies

* The HMF patient-db schema loaded (see patient-db project)
* Optional: set of recognised LINE element locations
* Optional: set of recognised fragile site locations

## Usage

```
java -jar sv-analysis-x.y-with-dependencies.jar
    -db_url "mysql://localhost/hmfpatients"
    -db_user username
    -db_pass password
    -sample SAMPLE_ID
``` 
