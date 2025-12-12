# Compar

A regression testing tool, comparing sample output across pipeline runs.

## Usage

```
java -jar compar.jar \
   -sample SAMPLE_T \
   -germline_sample SAMPLE_R \
   -categories ALL \
   -match_level REPORTABLE \
   -sample_dir_ref /path_to_sample_data/run_01/
   -sample_dir_new /path_to_sample_data/run_02/
   -output_dir /output_dir/ 
```

## Configuration
The key configuration values to set are:
- the sample(s) to compare
- the categories to compare - each of these will map to specific pipeline output files
- the source of data - either the MySQL hmf_patients DB or pipeline output files 
 
### Required configuration

| Argument                          | Description                                                                                                                                |
|-----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| sample                            | Tumor sample ID, OR                                                                                                                        |
| sample_id_file                    | File with column header SampleId and then list of sample IDs, optional Ref and New sample mappings and germline sample IDs (see examples)  |
| categories                        | 'ALL', 'PANEL', or otherwise specify a comma-separated list                                                                                |
| match_level                       | REPORTABLE (default) or DETAILED                                                                                                           |
| sample_data_ref & sample_data_new | Sample root directory for pipeline output                                                                                                  |
| TOOL_dir_ref & TOOL_dir_new **    | Tool path overrides - each pipeline tool directory eg 'linx_dir_ref' - relative path to 'sample_dir' if specified, otherwise absolute path |
| db_source_ref & db_source_new     | DB connection details for ref and new sample data - see format below                                                                       |
| output_dir                        | Path for output file                                                                                                                       |

** set of tools are: linx, linx_germline, purple, chord, cuppa, lilac, peach, virus (i.e. virus-interpreter), snp_genotype, tumor_flagstat, germline_flagstat, tumor_bam_metrics and germline_bam_metrics.

The available categories are: PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT, GERMLINE_DELETION, GERMLINE_SV, FUSION, DISRUPTION, CUPPA, 
CHORD, LILAC, PEACH, VIRUS, TUMOR_FLAGSTAT, GERMLINE_FLAGSTAT, TUMOR_BAM_METRICS, GERMLINE_BAM_METRICS, SNP_GENOTYPE, COPY_NUMBER, GENE_COPY_NUMBER,
CDR3_SEQUENCE, CDR3_LOCUS_SUMMARY, TELOMERE_LENGTH, V_CHORD.

The category PANEL is equivalent to PURITY, DRIVER, SOMATIC_VARIANT, FUSION, DISRUPTION, TUMOR_FLAGSTAT, TUMOR_BAM_METRICS and SNP_GENOTYPE, V_CHORD.


### Optional configuration

| Argument                                                | Description                                                                                                                                  |
|---------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| germline_sample                                         | Germline sample ID. Defaults to tumor sample ID with "-ref" appended                                                                         |
| output_id                                               | Outfile file suffix                                                                                                                          |
| driver_gene_panel                                       | Used to check alternate transcript changes and to limit analysis of somatics and gene copy number comparisons                                |
| restrict_to_drivers                                     | Limit analysis to genes within the panel                                                                                                     |
| write_detailed                                          | Write a file per compared category                                                                                                           |
| somatic_unfiltered_vcf_ref & somatic_unfiltered_vcf_new | VCF of unfiltered somatic variants (i.e. SAGE) for detecting filtering reason                                                                |
| liftover                                                | Apply liftover to relevant fields for pipeline run comparison across reference genome versions (V37/V38)                                     |
| include_matches                                         | Also include matching entries in output file(s)                                                                                              |
| pipeline_format_ref & pipeline_format_new               | Format for default tool directory derivation from sample directory. Default: OA_V2_3. Options: OA_V2_0, OA_V2_2, OA_V2_3, PIP5_V6_0, DB_V6_0 |
| pipeline_format_file_ref & pipeline_format_file_new     | Config file for default tool directory derivation from sample directory.                                                                     |


### Sample ID Mappings
If the same patient has different sample IDs for different runs and these are used for all filenames, then specify these mappings in the sample ID file, eg:
```
sample_id_mappings.csv
SampleId,RefSampleId,NewSampleId
COLO829T,COLO829_Ref,COLO829T_New
```
The same can be done for germline sample IDs.
```
sample_id_mappings.with_germline.csv
SampleId,GermlineSampleId,RefSampleId,RefGermlineSampleId,NewSampleId,NewGermlineSampleId
COLO829T,COLO829R,COLO829T_Ref,COLO829R_Ref,COLO829T_New,COLO829R_New
```

### File Sourced Data
Typically, set the 'sample_dir_ref' and 'sample_dir_new' to the REF and NEW sample root directories, which then contain each tool's 
output in a subdirectory as per the standard HMF pipeline.

Specify one or more tool directories to override the pipeline default paths.

| Category                                                     | Config Path Id    |
|--------------------------------------------------------------|-------------------|
| PURITY, SOMATIC_VARIANT, GERMLINE_VARIANT, GERMLINE_DELETION | purple_dir        |
| FUSION, DISRUPTION                                           | linx_dir          |
| GERMLINE_SV                                                  | linx_germline_dir |
| CUPPA                                                        | cuppa_dir         |
| CHORD                                                        | chord_dir         |
| LILAC                                                        | lilac_dir         |
| PEACH                                                        | peach_dir         |
| VIRUS                                                        | virus_dir         |
| CIDER                                                        | cider_dir         |
| TEAL                                                         | teal_dir          |
| V_CHORD                                                      | v_chord_dir       |

Wildcards '*' can be used in place of sampleIds, in which case Compar will replace the wildcard with the sampleId for each path.
Similarly, '$' can be used in place of germline sample IDs.

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

will load run 01 data Purple data from /path_to_sample_data/run_01/SAMPLE_ID/purple/.

It's possible to control the assumed pipeline output format for deriving the default tool paths from the sample data directories.
The `pipeline_format_ref` and/or `pipeline_format_new` arguments can be set to an older version of OncoAnalyser (e.g. `OA_V2_0`), our legacy pipeline5 WiGiTS implementation (`PIP5_V6_0`) or the format of the Hartwig Medical Database (`DB_V6_0`).
Alternatively, this format can be set in a config file such as [this](../hmf-common/src/test/resources/pipeline/completeToolDirectoryConfig.tsv) or [this](../hmf-common/src/test/resources/pipeline/partialToolDirectoryConfig.tsv)
by using the `pipeline_format_file_ref` and/or `pipeline_format_file_new` arguments.

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

| Field                         | Match Type & Thresholds |
|-------------------------------|-------------------------|
| qcStatus                      | Exact                   |
| gender                        | Exact                   |
| germlineAberration            | Exact                   |
| fitMethod                     | Exact                   |
| msStatus                      | Exact                   |
| tmbStatus                     | Exact                   |
| tmlStatus                     | Exact                   |
| purity                        | Threshold  [0.02]       |
| ploidy                        | Threshold  [0.1]        |
| contamination                 | Threshold  [0.005]      |
| tmbPerMb                      | Threshold  [0.1, 5%]    |
| msIndelsPerMb                 | Threshold  [0.1, 5%]    |
| tml                           | Threshold  [1, 5%]      |
| copyNumberSegments            | Threshold  [5, 20%]     |
| unsupportedCopyNumberSegments | Threshold  [5, 20%]     |
| svTmb                         | Threshold  [2, 5%]      |

### Somatic Variant
Data key: SampleId, Chromosome, Position, Ref, Alt and VariantType (SNP/MNP/INDEL/UNDEFINED)

| Field                      | Match Type & Thresholds |
|----------------------------|-------------------------|
| reported                   | Exact                   |
| filter                     | Exact                   |
| gene                       | Exact                   |
| canonicalEffect            | Exact                   |
| canonicalCodingEffect      | Exact                   |
| canonicalHgvsCodingImpact  | Exact                   |
| canonicalHgvsProteinImpact | Exact                   |
| otherTranscriptEffects     | Exact                   |
| tier                       | Exact                   |
| hotspot                    | Exact                   |
| biallelic                  | Exact                   |
| qual                       | Threshold [20, 20%]     |
| subclonalLikelihood        | Threshold [0.6]         |
| hasLPS                     | Exact                   |
| variantCopyNumber          | Threshold [0.3, 15%]    |
| tumorSupportingReadCount   | Threshold [1, 20%]      |
| tumorTotalReadCount        | Threshold [1, 20%]      |
| purityAdjustedVaf          | Threshold [0.2]         |

### Germline Variant
Data key: SampleId, Chromosome, Position, Ref, Alt and VariantType (SNP/MNP/INDEL/UNDEFINED)

| Field                      | Match Type & Thresholds |
|----------------------------|-------------------------|
| reported                   | Exact                   |
| filter                     | Exact                   |
| gene                       | Exact                   |
| canonicalEffect            | Exact                   |
| canonicalCodingEffect      | Exact                   |
| canonicalHgvsCodingImpact  | Exact                   |
| canonicalHgvsProteinImpact | Exact                   |
| otherTranscriptEffects     | Exact                   |
| tier                       | Exact                   |
| hotspot                    | Exact                   |
| biallelic                  | Exact                   |
| qual                       | Threshold [20, 20%]     |
| variantCopyNumber (tumor)  | Threshold [0.3, 15%]    |
| tumorSupportingReadCount   | Threshold [1, 20%]      |
| tumorTotalReadCount        | Threshold [1, 20%]      |
| purityAdjustedVaf          | Threshold [0.2]         |

### Germline Deletion
Data key: SampleId, Gene

| Field              | Match Type & Thresholds |
|--------------------|-------------------------|
| reported           | Exact                   |
| germlineStatus     | Exact                   |
| tumorStatus        | Exact                   |
| germlineCopyNumber | Threshold [0.2, 10%]    |
| tumorCopyNumber    | Threshold [0.2, 10%]    |
| chromosome         | Exact                   |
| chromosomeBand     | Exact                   |

### Drivers (Linx and Purple)
Data key: SampleId, GeneId, TranscriptId, Driver-type

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| likelihoodMethod | Exact                   |
| driverLikelihood | Threshold [0.1]         |
| minCopyNumber    | Threshold [0.3, 15%]    |
| maxCopyNumber    | Threshold [0.3, 15%]    |
| chromosome       | Exact                   |
| chromosomeBand   | Exact                   |

### Fusions
Data key: SampleId, Fusion name

| Field               | Match Type & Thresholds |
|---------------------|-------------------------|
| reported            | Exact                   |
| reportedType        | Exact                   |
| phased              | Exact                   |
| likelihood          | Exact                   |
| fusedExonUp         | Exact                   |
| fusedTranscriptUp   | Exact                   |
| fusedExonDown       | Exact                   |
| fusedTranscriptDown | Exact                   |
| chainLinks          | Exact                   |
| chainTerminated     | Exact                   |
| domainsKept         | Exact                   |
| domainsLost         | Exact                   |
| junctionCopyNumber  | Threshold [0.5, 20%]    |

### Disruptions
Data key: SampleId, Gene

Per breakend key: SV coordinates (chromosome, position, orientation), TranscriptId (if not canonical)

| Field              | Match Type & Thresholds |
|--------------------|-------------------------|
| reported           | Exact                   |
| regionType         | Exact                   |
| codingType         | Exact                   |
| nextSpliceExonRank | Exact                   |

### Germline SV
Data key: SampleId, Gene

Per breakend key: SV coordinates (chromosome, position, orientation), TranscriptId (if not canonical)

| Field              | Match Type & Thresholds |
|--------------------|-------------------------|
| reported           | Exact                   |
| regionType         | Exact                   |
| codingType         | Exact                   |
| nextSpliceExonRank | Exact                   |

### Cuppa
Data key: SampleId, ClassifierName, DataType

| Field         | Match Type & Thresholds |
|---------------|-------------------------|
| topCancerType | Exact                   |
| probability   | Threshold [0.1]         |

### Chord
Data key: SampleId

| Field    | Match Type & Thresholds |
|----------|-------------------------|
| BRCA1    | Threshold [0.1]         |
| BRCA2    | Threshold [0.1]         |
| status   | Exact                   |
| type     | Exact                   |
| hrdScore | Threshold [0.1]         |

### Lilac
Data key: SampleId

| Field                       | Match Type & Thresholds                                  |
|-----------------------------|----------------------------------------------------------|
| status                      | Exact                                                    |
| alleles                     | Exact - checks all 6 alleles match                       |
| somaticMissense             | Threshold [0.4, 10%] - checks per allele                 |
| somaticNonsenseOrFrameshift | Threshold [0.4, 10%] - checks per allele                 |
| somaticSplice               | Threshold [0.4, 10%] - checks per allele                 |
| somaticInframeIndel         | Threshold [0.4, 10%] - checks per allele                 |
| somaticSynonymous           | Threshold [0.4, 10%] - checks per allele (DETAILED only) |
| refTotalFragments           | Threshold [10, 1%] - checks per allele (DETAILED only)   |
| tumorTotalFragments         | Threshold [10, 1%] - checks per allele (DETAILED only)   |
| tumorCopyNumber             | Threshold [0.5, 15%] - checks per allele                 |
| totalFragments              | Threshold [10, 1%]                                       |
| fittedFragments             | Threshold [10, 1%]                                       |
| discardedAlignmentFragments | Threshold [10, 1%]                                       |
| discardedIndels             | Threshold [10, 1%]                                       |
| hlaYAllele                  | Exact                                                    |

### Peach
Data key: SampleId, Gene, HaplotypeName

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| alleleCount      | Exact                   |
| function         | Exact                   |
| drugs            | Exact                   | 
| prescriptionUrls | Exact                   | 

### Virus
Data key: SampleId, Name

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| reported         | Exact                   |
| integrations     | Threshold [20%]         |
| meanCoverage     | Threshold [15%]         |
| driverLikelihood | Exact                   |

### Flagstat
For both tumor and germline sample

Data key: SampleId

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| mappedProportion | Threshold [0.01]        |

### Tumor BAM Metrics
Data key: SampleId

| Field               | Match Type & Thresholds |
|---------------------|-------------------------|
| duplicatePercentage | Threshold [0.05]        |
| percentage30X       | Threshold [0.03]        |
| percentage60X       | Threshold [0.03]        |

### Germline BAM Metrics
Data key: SampleId

| Field               | Match Type & Thresholds |
|---------------------|-------------------------|
| duplicatePercentage | Threshold [0.05]        |
| percentage10X       | Threshold [0.03]        |
| percentage20X       | Threshold [0.03]        |

### SNP Genotype
Data key: SampleId, Chromosome, Position, Ref

| Field    | Match Type & Thresholds |
|----------|-------------------------|
| alt      | Exact                   |
| genotype | Exact                   |

### CDR3 Sequence (Cider)
Data key: SampleId, Cdr3AA, Cdr3Seq

| Field  | Match Type & Thresholds |
|--------|-------------------------|
| Filter | Exact                   |
| Locus  | Exact                   |

### CDR3 Locus Summary (Cider)
Data key: SampleId, Locus

| Field         | Match Type & Thresholds |
|---------------|-------------------------|
| passSequences | Threshold [5%]          | 

### Telomere Lengths
Data key: SampleId, Type (tumor or ref)

| Field          | Match Type & Thresholds |
|----------------|-------------------------|
| telomereLength | Threshold [5%]          | 

### vChord
Data key: SampleId

| Field                 | Match Type & Thresholds |
|-----------------------|-------------------------|
| breastCancerHrdScore  | Threshold [0.1]         | 
| ovarianCancerHrdScore | Threshold [0.1]         | 
| pancreaticCancerScore | Threshold [0.1]         | 
| prostateCancerScore   | Threshold [0.1]         | 
| otherCancerScore      | Threshold [0.1]         | 

### Copy Number
Only runs in DETAILED mode.

Data key: SampleId, Chromosome, StartPosition, EndPosition

| Field                | Match Type & Thresholds |
|----------------------|-------------------------|
| copyNumber           | Threshold [0.5, 15%]    |
| majorAlleleCopyNumer | Threshold [0.5, 15%]    |
| method               | Exact                   |

### Gene Copy Number
Only runs in DETAILED mode.

Data key: SampleId, Gene

| Field         | Match Type & Thresholds |
|---------------|-------------------------|
| minCopyNumber | Threshold [0.5, 15%]    |
| maxCopyNumer  | Threshold [0.5, 15%]    |

## Version History and Download Links
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/compar-v1.3)
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/compar-v1.2)
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/compar-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/compar-v1.0)
