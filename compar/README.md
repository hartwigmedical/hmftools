# Compar

A regression testing tool, comparing sample output across pipeline runs.

## Usage

```
java -jar compar.jar \
   -sample SAMPLE_T \
   -categories ALL \
   -match_level REPORTABLE \
   -file_sources "RUN_01;sample_dir=./run_01/,RUN_02;sample_dir=./run_02/" \
   -output_dir /output_dir/ 
```

## Configuration
 
The following filters are applied to all variants

Filter | Description
---|---
sample | Tumor sample ID or
sample_id_file | File with column header SampleId and then list of sample IDs
categories | 'ALL', otherwise specify a comma-separated list from PURITY, COPY_NUMBER, DRIVER, SOMATIC_VARIANT, LINX_DATA, FUSION, DISRUPTION
match_level | REPORTABLE, KEY_FIELDS or DETAILED
file_sources | List of sources and their file locations - see format below
db_sources |  List of sources and their DB locations - see format below
output_dir | Path for output file
output_id | Optional: outfile file suffix
source_sample_mappings | Optional: suffix to add to each sampleId by source, eg "run_01=_01,run_02=_02"

### File Sourced Data
Specify 'file_sources' config with a comma-separated list of the follow:
- SourceName
- optional sample directory - path to directory containing sample files
- optional Linx, Purple and Somatic file directories, with relative path used if sample directory is specified

Example 1
```
file_sources="RUN_01;sample_dir=/path_to_sample_data/run_01/,RUN_02;sample_dir=/path_to_sample_data/run_02/"
```

will load run 01 data from /path_to_sample_data/run_01/ and run 02 data from /path_to_sample_data/run_02/

Example 2
```
file_sources="RUN_01;sample_dir=/path_to_sample_data/run_01/;linx_dir=linx;purple_dir=purple,RUN_02 etc"
```

will load run 01 data Linx data from /path_to_sample_data/run_01/linx/ and Purple data from /path_to_sample_data/run_01/purple/ 

Example 3
```
file_sources="RUN_01;linx_dir=/path_to_sample_data/run_01/linx;purple_dir=/path_to_sample_data/run_01/purple/,RUN_02 etc"
```

will load run 01 data Linx data from /path_to_sample_data/run_01/linx/ and Purple data from /path_to_sample_data/run_01/purple/ 

### Database Sourced Data
Specify 'db_sources' config with a comma-separated list of the follow:
- SourceName;DbURL;DbUser;DbPassword

Example:
```
db_sources="PROD;mysql://localhost/prod;user1;pass1,TEST;mysql://localhost/test;user1;pass1"
```


## Data Categories
The current types for comparison are:
- Purple purity
- Purple copy number
- Purple somatic variants
- Purple structural variants
- Linx SV annotations and clustering
- Linx fusions
- Linx disruptions
- drivers (Linx and Purple)


## Match Levels
Comparison can be performed at 3 levels:
- REPORTABLE - only items with a reportable status are compared, and then only written as a mismatch if their reportable status diffs
- MODERATE - key fields are compared for each category - see below for details
- DETAILED - not yet implemented


## Version History and Download Links
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/compar-v1.0)
