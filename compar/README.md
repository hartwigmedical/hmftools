# Compar

A regression testing tool, comparing sample output across pipeline runs.

## Usage

```
java -jar compar.jar \
   -sample SAMPLE_T \
   -germline_sample SAMPLE_R \
   -categories ALL \
   -match_level REPORTABLE \
   -sample_dir_old /path_to_sample_data/run_01/
   -sample_dir_new /path_to_sample_data/run_02/
   -output_dir /output_dir/ 
```

## Configuration
The key configuration values to set are:
- the sample(s) to compare
- the categories to compare - each of these will map to specific pipeline output files
- the source of data - either the MySQL hmf_patients DB or pipeline output files 
 
### Required configuration

| Argument                          | Description                                                                                                                                          |
|-----------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| sample                            | Tumor sample ID, OR                                                                                                                                  |
| sample_id_file                    | File with column header SampleId and then list of sample IDs, optional Old and New sample mappings and refrence (germline) sample IDs (see examples) |
| categories                        | 'ALL', 'PANEL', or otherwise specify a comma-separated list                                                                                          |
| match_level                       | REPORTABLE (default) or DETAILED                                                                                                                     |
| sample_data_old & sample_data_new | Sample root directory for pipeline output                                                                                                            |
| TOOL_dir_old & TOOL_dir_new **    | Tool path overrides - each pipeline tool directory eg 'linx_dir_old' - relative path to 'sample_dir' if specified, otherwise absolute path           |
| db_source_old & db_source_new     | DB connection details for old and new sample data - see format below                                                                                 |
| output_dir                        | Path for output file                                                                                                                                 |

** set of tools are: linx, linx_germline, purple, chord, cuppa, isofox, lilac, peach, virus (i.e. virus-interpreter), sigs, snp_genotype, tumor_flagstat, germline_flagstat, tumor_bam_metrics and germline_bam_metrics.

The available categories are: PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT, GERMLINE_AMP_DEL, GERMLINE_SV, FUSION, DISRUPTION, CUPPA, CUPPA_IMAGE,
CHORD, LILAC_QC, LILAC_ALLELE, PEACH, VIRUS, TUMOR_FLAGSTAT, GERMLINE_FLAGSTAT, TUMOR_BAM_METRICS, GERMLINE_BAM_METRICS, SNP_GENOTYPE, COPY_NUMBER, GENE_COPY_NUMBER,
CDR3_SEQUENCE, CDR3_LOCUS_SUMMARY, TELOMERE_LENGTH, V_CHORD, SIGS, RNA_SUMMARY, RNA_GENE_DATA, RNA_TRANSCRIPT_DATA, NOVEL_SPLICE_JUNCTION, RNA_FUSION.

The category PANEL is equivalent to PURITY, DRIVER, SOMATIC_VARIANT, FUSION, DISRUPTION, TUMOR_FLAGSTAT, TUMOR_BAM_METRICS and SNP_GENOTYPE, V_CHORD.


### Optional configuration

| Argument                                                | Description                                                                                                                                     |
|---------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| germline_sample                                         | Germline sample ID. Defaults to tumor sample ID with "-ref" appended                                                                            |
| output_id                                               | Outfile file suffix                                                                                                                             |
| driver_gene_panel                                       | Used to check alternate transcript changes and to limit analysis of somatics and gene copy number comparisons                                   |
| restrict_to_drivers                                     | Limit analysis to genes within the panel                                                                                                        |
| write_detailed                                          | Write a file per compared category                                                                                                              |
| somatic_unfiltered_vcf_old & somatic_unfiltered_vcf_new | VCF of unfiltered somatic variants (i.e. SAGE) for detecting filtering reason                                                                   |
| liftover                                                | Apply liftover to relevant fields for pipeline run comparison across reference genome versions (V37/V38)                                        |
| include_matches                                         | Also include matching entries in output file(s)                                                                                                 |
| pipeline_format_old & pipeline_format_new               | Format for default tool directory derivation from sample directory. Default: OA_V2_3. Options: OA_V2_0, OA_V2_2, OA_V2_3, PIP5_V6_0, DB_V6_0    |
| pipeline_format_file_old & pipeline_format_file_new     | Config file for default tool directory derivation from sample directory.                                                                        |
| field_config_file                                       | Field config file (see [Field Config Files](#field-config-files)) used to override the default per-field `Compared`/threshold settings          |
| strict_field_config                                     | Requires `field_config_file` to explicitly set every compared field for every category being run; Compar exits with an error if any are missing |


### Sample ID Mappings
If the same patient has different sample IDs for different runs and these are used for all filenames, then specify these mappings in the sample ID file, eg:
```
sample_id_mappings.csv
SampleId,OldSampleId,NewSampleId
COLO829T,COLO829_Old,COLO829T_New
```
The same can be done for germline sample IDs.
```
sample_id_mappings.with_germline.csv
SampleId,ReferenceId,OldSampleId,OldReferenceId,NewSampleId,NewReferenceId
COLO829T,COLO829R,COLO829T_Old,COLO829R_Old,COLO829T_New,COLO829R_New
```

### File Sourced Data
Typically, set the 'sample_dir_old' and 'sample_dir_new' to the OLD and NEW sample root directories, which then contain each tool's 
output in a subdirectory as per the standard HMF pipeline.

Specify one or more tool directories to override the pipeline default paths.

| Category                                                                           | Config Path Id    |
|------------------------------------------------------------------------------------|-------------------|
| PURITY, SOMATIC_VARIANT, GERMLINE_VARIANT, GERMLINE_AMP_DEL                        | purple_dir        |
| FUSION, DISRUPTION                                                                 | linx_dir          |
| GERMLINE_SV                                                                        | linx_germline_dir |
| CUPPA                                                                              | cuppa_dir         |
| CHORD                                                                              | chord_dir         |
| LILAC_QC, LILAC_ALLELE                                                             | lilac_dir         |
| PEACH                                                                              | peach_dir         |
| VIRUS                                                                              | virus_dir         |
| CIDER                                                                              | cider_dir         |
| TEAL                                                                               | teal_dir          |
| V_CHORD                                                                            | v_chord_dir       |
| SIGS                                                                               | sigs_dir          |
| RNA_SUMMARY, RNA_GENE_DATA, RNA_TRANSCRIPT_DATA, NOVEL_SPLICE_JUNCTION, RNA_FUSION | isofox_dir        |

Wildcards '*' can be used in place of sampleIds, in which case Compar will replace the wildcard with the sampleId for each path.
Similarly, '$' can be used in place of germline sample IDs.

Example 1
```
purple_old="sample_dir=/path_to_sample_data/run_01/"
purple_new="sample_dir=/path_to_sample_data/run_02/"
```

will load old run 01 data from /path_to_sample_data/run_01/ and new run 02 data from /path_to_sample_data/run_02/

Example 2
```
"purple_old=/path_to_purple_data/run_01/*/purple/"
"purple_new=/path_to_purple_data/run_02/*/purple/"
```

will load run 01 data Purple data from /path_to_sample_data/run_01/SAMPLE_ID/purple/.

It's possible to control the assumed pipeline output format for deriving the default tool paths from the sample data directories.
The `pipeline_format_old` and/or `pipeline_format_new` arguments can be set to an older version of OncoAnalyser (e.g. `OA_V2_0`), our legacy pipeline5 WiGiTS implementation (`PIP5_V6_0`) or the format of the Hartwig Medical Database (`DB_V6_0`).
Alternatively, this format can be set in a config file such as [this](../hmf-common/src/test/resources/pipeline/completeToolDirectoryConfig.tsv) or [this](../hmf-common/src/test/resources/pipeline/partialToolDirectoryConfig.tsv)
by using the `pipeline_format_file_old` and/or `pipeline_format_file_new` arguments.

### Database Sourced Data
Specify 'db_sources' config with a comma-separated list of the follow:
- DbURL;DbUser;DbPassword

Example:
```
db_source_old="mysql://localhost/prod;user1;pass1"
db_source_new="mysql://localhost/test;user1;pass1"
```


## Data Categories, Fields and Thresholds
Each data type that is compared is described below. 
Differences in field values are considered one of the following ways
- an exact match, eg a string, boolean, string list or breakend value
- absolute difference vs a threshold 
- percentage difference vs a threshold
- absolute and percentage differences vs 2 thresholds, requiring both to be exceeded

### Field Config Files
The `Compared` setting and default thresholds for every field are defined per category in code, and are exported for reference to:
- [field.config.compar.tsv](src/main/resources/field_config/reportable/field.config.compar.tsv) for the default `-match_level REPORTABLE`
- [field.config.compar.tsv](src/main/resources/field_config/detailed/field.config.compar.tsv) for `-match_level DETAILED`, which also includes the extra categories and fields that are only compared in detailed mode

Each file lists one row per field, for every category that is run at that match level:

| Column            | Description                                                                                      |
|-------------------|--------------------------------------------------------------------------------------------------|
| Category          | The CategoryType the field belongs to                                                            |
| Field             | Field name, as it appears in mismatch output                                                     |
| FieldType         | string, stringList, boolean, int, long, double, breakend or pixel                                |
| Compared          | Whether the field is included in diff comparisons (some are context-only, or DETAILED mode only) |
| AbsoluteThreshold | Absolute difference threshold, or 'none'                                                         |
| PercentThreshold  | Percentage difference threshold, or 'none'                                                       |

These defaults can be overridden with `-field_config_file`, pointing to a TSV file in the same format containing only the rows to be overridden.
If `-strict_field_config` is also set, every compared field for every category being run must have an entry in the file, and Compar will exit
with an error if any are missing.

Every Compar run also writes out the field config it actually used (i.e. after applying any overrides) to `field.config.compar.tsv` in the
output directory. This is a record of exactly what was compared for that run, and can itself be used as a `-field_config_file` input to repeat
the comparison later with the same settings.

### Purity
Data key: SampleId

| Field                         | Match Type & Thresholds |
|-------------------------------|-------------------------|
| QcStatus                      | Exact                   |
| Gender                        | Exact                   |
| GermlineAberrations           | Exact                   |
| FitMethod                     | Exact                   |
| MsStatus                      | Exact                   |
| TmbStatus                     | Exact                   |
| TmlStatus                     | Exact                   |
| Purity                        | Threshold [0.04]        |
| Ploidy                        | Threshold [0.1]         |
| Contamination                 | Threshold [0.005]       |
| TmbPerMb                      | Threshold [0.1, 5%]     |
| MsIndelsPerMb                 | Threshold [0.1, 5%]     |
| Tml                           | Threshold [1, 5%]       |
| CopyNumberSegments            | Threshold [5, 20%]      |
| UnsupportedCopyNumberSegments | Threshold [5, 20%]      |
| SvTmb                         | Threshold [5, 5%]       |
| TincLevel                     | Threshold [0.1]         |

### Somatic Variant
Data key: SampleId, Chromosome, Position, Ref, Alt and VariantType (SNP/MNP/INDEL/UNDEFINED)

| Field                    | Match Type & Thresholds |
|--------------------------|-------------------------|
| Reported                 | Exact                   |
| filter                   | Exact                   |
| Gene                     | Exact                   |
| CanonicalEffect          | Exact                   |
| CanonicalCodingEffect    | Exact                   |
| CanonicalHgvsCoding      | Exact                   |
| CanonicalHgvsProtein     | Exact                   |
| OtherReportedEffects     | Exact                   |
| Tier                     | Exact                   |
| Hotspot                  | Exact                   |
| Biallelic                | Exact                   |
| BiallelicProb            | Threshold [0.3]         |
| Qual                     | Threshold [50, 20%]     |
| SubclonalLikelihood      | Threshold [0.6]         |
| HasLPS                   | Exact                   |
| VariantCopyNumber        | Threshold [0.3, 30%]    |
| TumorSupportingReadCount | Threshold [1, 20%]      |
| TumorTotalReadCount      | Threshold [1, 20%]      |
| PurityAdjustedVaf        | Threshold [0.2]         |

### Germline Variant
Data key: SampleId, Chromosome, Position, Ref, Alt and VariantType (SNP/MNP/INDEL/UNDEFINED)

| Field                    | Match Type & Thresholds |
|--------------------------|-------------------------|
| Reported                 | Exact                   |
| filter                   | Exact                   |
| Gene                     | Exact                   |
| CanonicalEffect          | Exact                   |
| CanonicalCodingEffect    | Exact                   |
| CanonicalHgvsCoding      | Exact                   |
| CanonicalHgvsProtein     | Exact                   |
| OtherReportedEffects     | Exact                   |
| Tier                     | Exact                   |
| Hotspot                  | Exact                   |
| Qual                     | Threshold [50, 20%]     |
| VariantCopyNumber        | Threshold [0.3, 30%]    |
| TumorSupportingReadCount | Threshold [1, 20%]      |
| TumorTotalReadCount      | Threshold [1, 20%]      |
| PurityAdjustedVaf        | Threshold [0.2]         |

### Germline Amp-Del
Data key: SampleId, Gene

| Field              | Match Type & Thresholds |
|--------------------|-------------------------|
| Reported           | Exact                   |
| GermlineStatus     | Exact                   |
| TumorStatus        | Exact                   |
| GermlineCopyNumber | Threshold [0.2, 10%]    |
| TumorCopyNumber    | Threshold [0.2, 10%]    |
| Chromosome         | Exact                   |
| ChromosomeBand     | Exact                   |

### Drivers (Linx and Purple)
Data key: SampleId, GeneId, TranscriptId, Driver-type

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| LikelihoodMethod | Exact                   |
| Likelihood       | Threshold [0.1]         |
| MinCopyNumber    | Threshold [0.3, 15%]    |
| MaxCopyNumber    | Threshold [0.3, 15%]    |
| Chromosome       | Exact                   |
| ChromosomeBand   | Exact                   |

### Fusions
Data key: SampleId, Fusion name

| Field           | Match Type & Thresholds |
|-----------------|-------------------------|
| Reported        | Exact                   |
| ReportedType    | Exact                   |
| Phased          | Exact                   |
| Likelihood      | Exact                   |
| FusedExonUp     | Exact                   |
| FusedExonDown   | Exact                   |
| ChainLinks      | Exact                   |
| ChainTerminated | Exact                   |
| DomainsKept     | Exact                   |
| DomainsLost     | Exact                   |

### Disruptions
Data key: SampleId, Gene

| Field    | Match Type & Thresholds |
|----------|-------------------------|
| Breakend | Exact - see below       |

A diff is raised per gene if, for any matched breakend the reported status, region type, coding type, or next splice exon rank differs between the two runs,
or a breakend/SV present in one run has no match in the other.

### Germline SV
Data key: SampleId, Gene

| Field    | Match Type & Thresholds |
|----------|-------------------------|
| Breakend | Exact - see below       |

Compared the same way as [Disruptions](#disruptions) above.

### Cuppa
Data key: SampleId, ClassifierName, DataType

| Field           | Match Type & Thresholds |
|-----------------|-------------------------|
| top_cancer_type | Exact                   |
| probability     | Threshold [0.1]         |

### Cuppa Image
Data key: SampleId, Dimensions

| Field      | Match Type & Thresholds                                                    |
|------------|----------------------------------------------------------------------------|
| Dimensions | Exact                                                                      |
| Pixels     | Threshold [0%] - percentage of differing pixels between the two chart PNGs |

### Chord
Data key: SampleId

| Field  | Match Type & Thresholds |
|--------|-------------------------|
| BRCA1  | Threshold [0.1]         |
| BRCA2  | Threshold [0.1]         |
| Status | Exact                   |
| Type   | Exact                   |
| Score  | Threshold [0.1]         |

### Lilac QC
Data key: SampleId

| Field                       | Match Type & Thresholds            |
|-----------------------------|------------------------------------|
| Status                      | Exact                              |
| Alleles                     | Exact - checks all 6 alleles match |
| TotalFragments              | Threshold [10, 1%]                 |
| FittedFragments             | Threshold [10, 1%]                 |
| DiscardedAlignmentFragments | Threshold [10, 1%]                 |
| DiscardedIndels             | Threshold [10, 1%]                 |
| HlaYAllele                  | Exact                              |

### Lilac Allele
Data key: SampleId, Gene, Allele, occurrence index (starts at 1, to distinguish an allele called more than once for the same gene)

| Field                       | Match Type & Thresholds              |
|-----------------------------|--------------------------------------|
| SomaticMissense             | Threshold [0.4, 10%]                 |
| SomaticNonsenseOrFrameshift | Threshold [0.4, 10%]                 |
| SomaticSplice               | Threshold [0.4, 10%]                 |
| SomaticInframeIndel         | Threshold [0.4, 10%]                 |
| SomaticSynonymous           | Threshold [0.4, 10%] (DETAILED only) |
| RefTotal                    | Threshold [10, 1%] (DETAILED only)   |
| TumorTotal                  | Threshold [10, 1%] (DETAILED only)   |
| TumorCopyNumber             | Threshold [0.5, 15%]                 |

### Peach
Data key: SampleId, Gene, HaplotypeName

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| AlleleCount      | Exact                   |
| Function         | Exact                   |
| Drugs            | Exact                   |
| PrescriptionUrls | Exact                   |

### Virus
Data key: SampleId, Name

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| Reported         | Exact                   |
| Integrations     | Threshold [20%]         |
| MeanCoverage     | Threshold [15%]         |
| DriverLikelihood | Exact                   |

### Flagstat
For both tumor and germline sample

Data key: SampleId

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| MappedProportion | Threshold [0.01]        |

### Tumor BAM Metrics
Data key: SampleId

| Field               | Match Type & Thresholds |
|---------------------|-------------------------|
| DuplicatePercentage | Threshold [0.05]        |
| Percentage30X       | Threshold [0.03]        |
| Percentage60X       | Threshold [0.03]        |

### Germline BAM Metrics
Data key: SampleId

| Field               | Match Type & Thresholds |
|---------------------|-------------------------|
| DuplicatePercentage | Threshold [0.05]        |
| Percentage10X       | Threshold [0.03]        |
| Percentage20X       | Threshold [0.03]        |

### SNP Genotype
Data key: SampleId, Chromosome, Position, Ref

| Field       | Match Type & Thresholds |
|-------------|-------------------------|
| Alt         | Exact                   |
| Genotype    | Exact                   |
| VcfSampleId | Exact                   |

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
| PassSequences | Threshold [5%]          |

### Telomere Lengths
Data key: SampleId, Type (tumor or ref)

| Field          | Match Type & Thresholds |
|----------------|-------------------------|
| TelomereLength | Threshold [5%]          |

### vChord
Data key: SampleId

| Field                 | Match Type & Thresholds |
|-----------------------|-------------------------|
| BreastCancerHrdScore  | Threshold [0.1]         |
| OvarianCancerHrdScore | Threshold [0.1]         |
| PancreaticCancerScore | Threshold [0.1]         |
| ProstateCancerScore   | Threshold [0.1]         |
| OtherCancerScore      | Threshold [0.1]         |

### Sigs
Only runs in DETAILED mode.

Data key: SampleId, Signature

| Field   | Match Type & Thresholds |
|---------|-------------------------|
| Percent | Threshold [0.05]        |

### RNA Summary
Only runs in DETAILED mode.

Data key: SampleId

| Field                 | Match Type & Thresholds |
|-----------------------|-------------------------|
| QcStatus              | Exact                   |
| TotalFragments        | Threshold [10, 1%]      |
| DuplicateFragments    | Threshold [10, 1%]      |
| SplicedFragmentPerc   | Threshold [0.01, 5%]    |
| UnsplicedFragmentPerc | Threshold [0.01, 5%]    |
| AltFragmentPerc       | Threshold [0.01, 5%]    |
| ChimericFragmentPerc  | Threshold [0.01, 5%]    |
| SplicedGeneCount      | Threshold [10, 1%]      |
| ReadLength            | Exact                   |
| FragLength5th         | Threshold [5%]          |
| FragLength50th        | Threshold [5%]          |
| FragLength95th        | Threshold [5%]          |
| EnrichedGenePercent   | Threshold [0.01]        |
| MedianGCRatio         | Threshold [0.01]        |
| ForwardStrandPercent  | Threshold [0.01]        |

### RNA Gene Data
Only runs in DETAILED mode.

Data key: SampleId, GeneName

| Field              | Match Type & Thresholds |
|--------------------|-------------------------|
| SplicedFragments   | Threshold [10, 5%]      |
| UnsplicedFragments | Threshold [10, 5%]      |
| AdjTPM             | Threshold [5%]          |

### RNA Transcript Data
Only runs in DETAILED mode.

Data key: SampleId, TranscriptName

| Field    | Match Type & Thresholds |
|----------|-------------------------|
| GeneName | Exact                   |
| AdjTPM   | Threshold [5%]          |

### Novel Splice Junction
Only runs in DETAILED mode.

Data key: SampleId, GeneName, Chromosome, JunctionStart, JunctionEnd

| Field       | Match Type & Thresholds |
|-------------|-------------------------|
| Type        | Exact                   |
| FragCount   | Threshold [5, 5%]       |
| RegionStart | Exact                   |
| RegionEnd   | Exact                   |

### RNA Fusion
Only runs in DETAILED mode.

Data key: SampleId, FusionName, ChromosomeUp, PositionUp, ChromosomeDown, PositionDown

| Field           | Match Type & Thresholds |
|-----------------|-------------------------|
| KnownFusionType | Exact                   |
| JuncTypeUp      | Exact                   |
| JuncTypeDown    | Exact                   |
| SplitFrags      | Threshold [5, 5%]       |

### Copy Number
Only runs in DETAILED mode.

Data key: SampleId, Chromosome, StartPosition, EndPosition

| Field                 | Match Type & Thresholds |
|-----------------------|-------------------------|
| CopyNumber            | Threshold [0.5, 15%]    |
| MajorAlleleCopyNumber | Threshold [0.5, 15%]    |
| Method                | Exact                   |

### Gene Copy Number
Only runs in DETAILED mode.

Data key: SampleId, Gene

| Field         | Match Type & Thresholds |
|---------------|-------------------------|
| MinCopyNumber | Threshold [0.5, 15%]    |
| MaxCopyNumber | Threshold [0.5, 15%]    |
