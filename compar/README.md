# Compar

A regression testing tool, comparing sample output across pipeline runs.

## Usage

```
java -jar compar.jar \
   -sample SAMPLE_T \
   -categories ALL \
   -match_level REPORTABLE \
   -sample_dir_ref /path_to_sample_data/run_01/
   -sample_dir_new /path_to_sample_data/run_02/
   -output_dir /output_dir/ 
```

## Configuration
The key configuration values to set are:
- the sample(s) or to compare
- the categories to compare - each of these will map to specific pipeline output files
- the source of data - either the MySQL hmf_patients DB or pipeine output files 
 
### Required configuration

Filter | Description
---|---
sample | Tumor sample ID, OR
sample_id_file | File with column header SampleId and then list of sample IDs, optional Ref and New sample mappings (see example)
categories | 'ALL', otherwise specify a comma-separated list from PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT, GERMLINE_DELETION, GERMLINE_SV, FUSION, DISRUPTION, CUPPA, CHORD, LILAC
match_level | REPORTABLE (default) or DETAILED
sample_data_ref & sample_data_new | Sample root directory for pipeline output
TOOL_dir_ref & TOOL_dir_new ** | Tool path overrides - each pipeline tool directory eg 'linx_dir_ref' - relative path to 'sample_dir' if specified, otherwise absolute path
db_source_ref & db_source_new |  DB connection details for ref and new sample data - see format below
output_dir | Path for output file

** set of tools are: linx, linx_germline, purple, chord, cuppa and lilac

### Optional configuration

Filter | Description
---|---
output_id | Optional: outfile file suffix
driver_gene_panel | Used to check alternate transcript changes and to limit analysis of somatics and gene copy number comparisons
restrict_to_drivers | Limit analysis to genes within the panel
write_detailed | Write a file per compared category

### Sample ID Mappings
If the same patient has different sample IDs for different runs and these are used for all filenames, then specify these mappings in the sample ID file, eg:
```
sample_id_mappings.csv
SampleId,RefSampleId,NewSampleId
COLO829T,COLO829_Ref,COLO829T_New
```

### File Sourced Data
Typically set the 'sample_dir_ref' and 'sample_dir_new' to the REF and NEW sample root directories, which then contain each tool's 
output in a sub-directory as per the standard HMF pipeline.

Specify one or more tool directories to override the pipeline default paths.

Category | Config Path Id
---|---
PURITY, SOMATIC_VARIANT, GERMLINE_VARIANT, GERMLINE_DELETION | purple_dir
FUSION, DISRUPTION | linx_dir
GERMLINE_SV | linx_germline_dir
PURITY | purple_dir
CUPPA | cuppa_dir
CHORD | chord_dir
LILAC | lilac_dir

Wildcards '*' can be used in place of sampleIds, in which case Compar will replace the wildcard with the sampleId for each path.

Example 1
```
purple_ref="sample_dir=/path_to_sample_data/run_01/"
purple_new="sample_dir=/path_to_sample_data/run_02/"
```

will load reference run 01 data from /path_to_sample_data/run_01/ and new run 02 data from /path_to_sample_data/run_02/

Example 2
```
"purple_ref=/path_to_purple_data/run_01/*/purple/"
"purple_new=/path_to_purple_data/run_02/*/purple/"
```

will load run 01 data Purple data from /path_to_sample_data/run_01/SAMPLE_ID/purple/ 


### Database Sourced Data
Specify 'db_sources' config with a comma-separated list of the follow:
- DbURL;DbUser;DbPassword

Example:
```
db_source_ref="mysql://localhost/prod;user1;pass1"
db_source_new="mysql://localhost/test;user1;pass1"
```


## Data Categories, Fields and Thresholds
Each data type that is compared is described below. 
Differences in field values are considered one of the following ways
- an exact match, eg a string value or type
- absolute difference vs a threshold 
- percentage difference vs a threshold
- absolute and percentage differences vs 2 thresholds, requiring both to be exceeded
 

### Purity
Data key: SampleId

Field | Match Type & Thresholds 
---|---
qcStatus | Exact
gender | Exact
germlineAberration | Exact
fitMethod | Exact
msStatus | Exact
tmbStatus | Exact
tmlStatus | Exact
purity | Threshold  [0.02]
ploidy | Threshold  [0.1]
contamination | Threshold  [0.005]
tmbPerMb | Threshold  [0.1, 5%]
msIndelsPerMb | Threshold  [0.1, 5%]
tml | Threshold  [1, 5%]
copyNumberSegments | Threshold  [5, 20%]
unsupportedCopyNumberSegments | Threshold  [5, 20%]
svTmb | Threshold  [2, 5%]

### Somatic Variant
Data key: SampleId, Chromosome, Position, Ref and Alt

Field | Match Type & Thresholds 
---|---
reported | Exact
filter | Exact
gene | Exact
canonicalEffect | Exact
canonicalCodingEffect | Exact
canonicalHgvsCodingImpact | Exact
canonicalHgvsProteinImpact | Exact
otherTranscriptEffects | Exact
tier | Exact
hotspot | Exact
biallelic | Exact
qual | Threshold [20, 20%]
subclonalLikelihood | [0.6]
has LPS | true / false

### Germline Variant
Data key: SampleId, Chromosome, Position, Ref and Alt

Field | Match Type & Thresholds 
---|---
reported | Exact
filter | Exact
gene | Exact
canonicalEffect | Exact
canonicalCodingEffect | Exact
canonicalHgvsCodingImpact | Exact
canonicalHgvsProteinImpact | Exact
otherTranscriptEffects | Exact
tier | Exact
hotspot | Exact
biallelic | Exact
pathogenicity | Exact
pathogenic | Exact
qual | Threshold [20, 20%]

### Germline Deletion
Data key: SampleId, Gene

Field | Match Type & Thresholds 
---|---
reported | Exact
germlineStatus | Exact
tumorStatus | Exact
germlineCopyNumber | Threshold [0.2, 10%]
tumorCopyNumber | Threshold [0.2, 10%]

### Drivers (Linx and Purple)
Data key: SampleId, GeneId, TranscriptId, Driver-type

Field | Match Type & Thresholds 
---|---
likelihoodMethod | Exact
driverLikelihood | [0.1]
minCopyNumber | max(0.2, 10%) for AMPs and DELs - not currently checked

### Fusions
Data key: SampleId, Fusion name

Field | Match Type & Thresholds 
---|---
reported | Exact
reportedType | Exact
phased | Exact
likelihood | Exact
fusedExonUp | Exact
fusedExonDown | Exact
chainLinks | Exact
chainTerminated | Exact
domainsKept | Exact
domainsLost | Exact

### Disruptions
Data key: SampleId, SV coordinates (chromosome, position, orientation), Gene, TranscriptId

Field | Match Type & Thresholds 
---|---
reported | Exact
geneOrientation | Exact
regionType | Exact
codingContext | Exact
nextSpliceExonRank | Exact
undisruptedCopyNumber | max(0.2, 10%) - not currently checked

### Germline SV
Data key: SampleId, SV coordinates (chromosome, position, orientation), Gene

Field | Match Type & Thresholds 
---|---
reported | Exact
qual | Threshold [20, 20%]
germlineFragments | Threshold [5, 10%]

### Cuppa
Data key: SampleId, ClassifierName

Field | Match Type & Thresholds 
---|---
topRefCancerType | Exact
topRefValue | Threshold [10%]

### Chord
Data key: SampleId

Field | Match Type & Thresholds 
---|---
BRCA1 | Threshold [0.1]
BRCA2 | Threshold [0.1]
status | Exact
type | Exact
hrdScore | Threshold [0.1]

### Lilac
Data key: SampleId

Field | Match Type & Thresholds 
---|---
status | Exact
alleles | Exact - checks all 6 alleles match
somaticVariants | Exact - number of annotated somatic variants match

## Version History and Download Links
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/compar-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/compar-v1.0)
