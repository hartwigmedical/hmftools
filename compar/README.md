# Compar

A regression testing tool, comparing sample output across pipeline runs.

## Usage

```
java -jar compar.jar \
   -sample SAMPLE_T \
   -germline_sample SAMPLE_R \
   -categories ALL \
   -match_level REPORTABLE \
   -sample_dir_old /path_to_sample_data/run_01/ \
   -sample_dir_new /path_to_sample_data/run_02/ \
   -output_dir /output_dir/ 
```

## Configuration
The key configuration values to set are:
- the sample(s) to compare
- the categories to compare - each of these will map to specific pipeline output files
- the source of data - pipeline output files, or expected findings files in place of the old source
 
### Required configuration

| Argument                        | Description                                                                                                                                          |
|---------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| sample                          | Tumor sample ID, OR                                                                                                                                  |
| sample_id_file                  | File with column header SampleId and then list of sample IDs, optional Old and New sample mappings and refrence (germline) sample IDs (see examples) |
| categories                      | 'ALL', 'PANEL', or otherwise specify a comma-separated list                                                                                          |
| match_level                     | REPORTABLE (default) or DETAILED                                                                                                                     |
| sample_dir_old & sample_dir_new | Sample root directory for pipeline output                                                                                                            |
| TOOL_dir_old & TOOL_dir_new **  | Tool path overrides - each pipeline tool directory eg 'linx_dir_old' - relative path to 'sample_dir' if specified, otherwise absolute path           |
| output_dir                      | Path for output file                                                                                                                                 |

** set of tools are: linx, linx_germline, purple, chord, cuppa, isofox, lilac, peach, virus (i.e. virus-interpreter), sigs, snp_genotype, tumor_flagstat, germline_flagstat, tumor_bam_metrics and germline_bam_metrics.

The available categories are: PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT, GERMLINE_AMP_DEL, GERMLINE_SV, FUSION, DISRUPTION, CUPPA, CUPPA_IMAGE,
CHORD, LILAC_QC, LILAC_ALLELE, PEACH, VIRUS, TUMOR_FLAGSTAT, GERMLINE_FLAGSTAT, TUMOR_BAM_METRICS, GERMLINE_BAM_METRICS, SNP_GENOTYPE, COPY_NUMBER, GENE_COPY_NUMBER,
CDR3_SEQUENCE, CDR3_LOCUS_SUMMARY, TELOMERE_LENGTH, V_CHORD, SIGS, RNA_SUMMARY, RNA_GENE_DATA, RNA_TRANSCRIPT_DATA, RNA_NOVEL_SPLICE_JUNCTION, RNA_FUSION.

The category PANEL is equivalent to PURITY, DRIVER, SOMATIC_VARIANT, FUSION, DISRUPTION, TUMOR_FLAGSTAT, TUMOR_BAM_METRICS and SNP_GENOTYPE, V_CHORD.


### Optional configuration

| Argument                                                | Description                                                                                                                                           |
|---------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| germline_sample                                         | Germline sample ID. Defaults to tumor sample ID with "-ref" appended                                                                                  |
| output_id                                               | Outfile file suffix                                                                                                                                   |
| driver_gene_panel                                       | Used to check alternate transcript changes and to limit analysis of somatics and gene copy number comparisons                                         |
| restrict_to_drivers                                     | Limit analysis to genes within the panel                                                                                                              |
| write_detailed                                          | Write a file per compared category                                                                                                                    |
| somatic_unfiltered_vcf_old & somatic_unfiltered_vcf_new | VCF of unfiltered somatic variants (i.e. SAGE) for detecting filtering reason                                                                         |
| liftover                                                | Apply liftover to relevant fields for pipeline run comparison across reference genome versions (V37/V38)                                              |
| include_matches                                         | Also include matching entries in output file(s)                                                                                                       |
| pipeline_format_old & pipeline_format_new               | Format for default tool directory derivation from sample directory. Default: OA_V3_0. Options: OA_V2_0, OA_V2_2, OA_V2_3, OA_V3_0, PIP5_V6_0, DB_V6_0 |
| pipeline_format_file_old & pipeline_format_file_new     | Config file for default tool directory derivation from sample directory.                                                                              |
| field_config_file                                       | Field config file (see [Field Config Files](#field-config-files)) used to override the default per-field `Compared`/threshold settings                |
| strict_field_config                                     | Requires `field_config_file` to explicitly set every compared field for every category being run; Compar exits with an error if any are missing       |
| source_label_old & source_label_new                     | Labels for the old and new sources, written to the combined output (see [Curation Cache](#curation-cache))                                            |
| curation_cache                                          | Single curation TSV for the pipeline, used to acknowledge known differences (see [Curation Cache](#curation-cache))                                   |


### Sample ID Mappings
If the same patient has different sample IDs for different runs and these are used for all filenames, then specify these mappings in the sample ID file, eg:
```
sample_id_mappings.csv
SampleId,OldSampleId,NewSampleId
COLO829T,COLO829_Old,COLO829T_New
```
`SampleId` is written to the curation file and history is looked up by it. `OldSampleId` and `NewSampleId` only identify the source
files. See [Curation Cache](#curation-cache).

The same can be done for germline sample IDs.
```
sample_id_mappings.with_germline.csv
SampleId,ReferenceId,OldSampleId,OldReferenceId,NewSampleId,NewReferenceId
COLO829T,COLO829R,COLO829T_Old,COLO829R_Old,COLO829T_New,COLO829R_New
```

### Per-sample directories
If the sample root directories don't follow a common pattern, they can be specified per sample in the sample ID file with the optional
columns `OldSampleDir` and `NewSampleDir`. These override `sample_dir_old` and `sample_dir_new` respectively for that sample, while
any tool directory overrides and the pipeline format still apply. Leave a value empty to use the default sample directory.
```
sample_ids_with_dirs.csv
SampleId,OldSampleDir,NewSampleDir
COLO829T,/path_to_sample_data/run_01/COLO829T/,/path_to_sample_data/run_02/COLO829T/
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
| RNA_SUMMARY, RNA_GENE_DATA, RNA_TRANSCRIPT_DATA, RNA_NOVEL_SPLICE_JUNCTION, RNA_FUSION | isofox_dir        |

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

> [!NOTE]
> Compar supports output comparison from OncoAnalyser v3.0 onwards.

### Output Files
Output files are named the same for one sample or many: `compar.tsv` for all categories, or one `compar.<category>.tsv`
per category with `write_detailed`. When set, `output_id` follows `compar`, e.g.
`compar.compar-2026-10-01T10-35-26.somatic_variant.tsv`. Every row starts with its `SampleId`.

### Expected Findings Files
Set `expected_findings_files_old` in place of the old sample directory to compare a run with expected values. Separate files with ';'.

Columns:

| Column       | Use      | Description                                                                                |
|--------------|----------|--------------------------------------------------------------------------------------------|
| SampleId     | Optional | If no row has one, every row applies to every sample                                       |
| Category     | Required | A category from the table below                                                            |
| Key          | Required | The item's key, in a form from the table below                                             |
| Status       | Optional | `Present` or `Absent`. Left blank, it defaults to `Present`. Any other value stops the run |
| FieldValues  | Optional | `Field=Value` pairs separated by ';', empty for none                                       |
| *Field name* | Optional | One column per field, e.g. `MinCopyNumber`                                                 |
| Remark       | Optional | Free text, copied to the output and not compared                                           |

> [!NOTE]
> Give the values to compare either in a `FieldValues` column, as `Field=Value` pairs separated by ';', so one file can hold every
> category, or in one column per field, with one file per category.

Use `N/A` in a field column for a value that should not be compared.

| Category                  | Key                                                                   |
|---------------------------|-----------------------------------------------------------------------|
| PURITY                    | `PURITY`                                                              |
| SOMATIC_VARIANT           | `Chromosome:Position:Ref:Alt`                                         |
| SOMATIC_VARIANT           | `GeneName:HgvsProtein`, e.g. `BRAF:p.V600E`                           |
| SOMATIC_VARIANT           | `GeneName`, for an `Absent` row only                                  |
| GENE_COPY_NUMBER          | `GeneName`, e.g. `ERBB2`                                              |
| GENE_COPY_NUMBER          | `GeneName:TranscriptId`, e.g. `CDKN2A:ENST00000579755`                |
| RNA_NOVEL_SPLICE_JUNCTION | `GeneName:Chromosome:JunctionStart:JunctionEnd`                       |
| RNA_FUSION                | `GeneUp:GeneDown:ChromosomeUp:PositionUp:ChromosomeDown:PositionDown` |
| RNA_FUSION                | `GeneUp:GeneDown`                                                     |

- Field names are in [Data Categories, Fields and Thresholds](#data-categories-fields-and-thresholds).
- Every row counts as passing, and a row with no matching pipeline row is `OLD_ONLY` at REPORTABLE too. To compare reported status,
  give `Reported=true` (SOMATIC_VARIANT only).

#### Absent rows
An `Absent` row means the run must not find that item. Its field values are not compared.

- Only SOMATIC_VARIANT and RNA_FUSION take `Absent` rows, and the key can't have coordinates. Breaking either rule stops the run.
- A SOMATIC_VARIANT gene key, e.g. `KRAS`, matches every pipeline row in that gene. Any other key matches what a `Present` row with
  that key would.
- Only pipeline rows that no `Present` row matches are checked: reported rows at REPORTABLE, passing rows at DETAILED.

| When                            | Output row                                         |
|---------------------------------|----------------------------------------------------|
| A checked pipeline row matches  | `ABSENT_FOUND`, keyed by the pipeline row          |
| No checked pipeline row matches | `ABSENT_CONFIRMED`, even without `include_matches` |

#### By category

##### Purity
There is one record per sample. A sample with no expected findings row is not compared, and a sample with expected findings but no
pipeline output gives `INVALID_NEW`.

```tsv
SampleId	Category	Key	FieldValues	Remark
DEMO_001	PURITY	PURITY	Purity=0.65;Ploidy=3.1;MsStatus=MSS	
```

##### Somatic variants
- `GeneName:HgvsProtein` matches a call in that gene with the same canonical protein change, independent of the reference genome.
  PAVE's notation and single-letter amino acids are both accepted. A variant with no protein change, such as `p.?`, needs a genomic key.
- A variant the run filtered is `OLD_ONLY`, with a `Filter(PASS/PON)` difference when Compar finds its VCF record.

```tsv
SampleId	Category	Key	Status	FieldValues	Remark
DEMO_001	SOMATIC_VARIANT	chr7:140753336:A:T		AF=0.30	
DEMO_001	SOMATIC_VARIANT	EGFR:p.L858R		AF=0.25	
DEMO_001	SOMATIC_VARIANT	KRAS	Absent		wild type by ddPCR
```

##### Gene copy numbers
The key is the gene name. `MinCopyNumber`, `MaxCopyNumber` and `DriverType` are compared. `DriverType` is Purple's call, such as `AMP`,
`DEL` or `HET_DEL`, and `UNKNOWN` when it made none. The run only loads `driver_gene_panel` genes, so a listed gene outside the panel is
`OLD_ONLY`, and unlisted panel genes are not compared.

```tsv
SampleId	Category	Key	MinCopyNumber	MaxCopyNumber	DriverType	Remark
DEMO_001	GENE_COPY_NUMBER	ERBB2	N/A	N/A	AMP	amplified by FISH
DEMO_001	GENE_COPY_NUMBER	CDKN2A	0	N/A	DEL	
```

> [!NOTE]
> Providing a transcript ID with the gene name, as `GeneName:TranscriptId` (e.g. `CDKN2A:ENST00000579755`), compares that
> transcript's copy number only. Use one kind of key per gene.

##### RNA novel splice junctions
- The key is `GeneName:Chromosome:JunctionStart:JunctionEnd`, in the run's reference genome. A pipeline row must match all four.
- `Type` is the junction type as Isofox writes it, such as `SKIPPED_EXONS` or `NOVEL_EXON`. See
  [Isofox novel splice junctions](../isofox/README.md#7-counting-and-characterisation-of-novel-splice-junctions) for what each means.
- Every junction in the run that the file doesn't list is `NEW_ONLY`.

```tsv
SampleId	Category	Key	Type	FragmentCount	Remark
DEMO_001	RNA_NOVEL_SPLICE_JUNCTION	MET:chr7:1000:2000	SKIPPED_EXONS	40	confirmed by RT-PCR
DEMO_001	RNA_NOVEL_SPLICE_JUNCTION	TP53:chr17:3000:4000	NOVEL_EXON	N/A	
```

##### RNA fusions
- `GeneUp:GeneDown` needs any passing fusion between those genes, and its field values are not compared. If none passed, the row shows
  the filter of the unfiltered fusion with the most split fragments.
- An Absent `GeneUp:GeneDown` key counts every passing fusion between those genes, one `ABSENT_FOUND` row per breakend.
- A fusion with both a gene-pair key and a coordinate key is not compared, and RNA_FUSION is `INVALID_ERROR` for that sample.

```tsv
SampleId	Category	Key	Status	FieldValues	Remark
DEMO_001	RNA_FUSION	EML4:ALK:chr2:3000:chr2:4000		KnownType=KNOWN_PAIR;SplitFrags=12	confirmed by RT-PCR
DEMO_002	RNA_FUSION	TMPRSS2:ERG			
DEMO_002	RNA_FUSION	BCR:ABL1	Absent		negative by FISH
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
- [example.field_config.tsv](src/main/resources/field_config/reportable/example.field_config.tsv) for the default `-match_level REPORTABLE`
- [example.field_config.tsv](src/main/resources/field_config/detailed/example.field_config.tsv) for `-match_level DETAILED`, which also includes the extra categories and fields that are only compared in detailed mode

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

Optionally, a Compar run also writes out the field config it actually used (i.e. after applying any overrides) to a `.field_config.tsv` file in the
output directory. This is a record of exactly what was compared for that run, and can itself be used as a `-field_config_file` input to repeat
the comparison later with the same settings.

### Curation Cache

The curation cache is a TSV of reviewed differences, one per pipeline. Compar reads it to mark matching output rows and never writes
to it.

```bash
-source_label_old v34 -source_label_new v35
```

#### Source labels

Source labels identify each OLD/NEW comparison in the curation cache and appear in previous comparison details.
A label is semicolon-separated `key=value` fields in this order: `kind`, `version`, optional `reference`, optional `config`, e.g.
`kind=pipeline;version=v34;reference=GRCh38;config=standard`. A label in any other form is kept as plain text.
Labels must be non-blank, with no tabs or line breaks. Source labels and `-curation_cache` require `GENERIC` output.

#### Cache columns

A cache row is a row of the combined output, copied up to `NewSource`, followed by the curation columns. Compar needs these columns:
`SampleId`, `Category`, `MismatchType`, `Key`, `Differences`, `OldSource`, `NewSource`, `CurationType`, `Severity`, `Verdict`,
`CurationComment`, `CuratedBy` and `CuratedDate`. A missing column stops the run. `CurationStatus`, `OldValues` and `NewValues` may be
present and are ignored. `CurationComment` is free text and may be blank. A cache row applies to output rows with the same `SampleId`,
`Category` and `Key`.

#### Severity

`Severity` is the review label in each cache row. Allowed values are `NONE`, `CRITICAL`, `MEDIUM` and `LOW`. Any other value,
including a blank cell, stops the run.

#### Curation type

`CurationType` says how a cache row was made: `MANUAL` or `BATCH_INFERRED`. Any other value, including a blank cell, stops the run.

#### Verdict

`Verdict` is the reviewer's call on which side is right: `NEW_CORRECT`, `OLD_CORRECT`, `BOTH_CORRECT`, `MIXED`, `UNCLEAR` or `NONE`.
A blank cell means `NONE`. Any other value stops the run.

#### CuratedBy and CuratedDate

`CuratedBy` and `CuratedDate` record who last changed the cache row and when. Both are required on every row. `CuratedDate` is an ISO
date (`2026-10-20`) or an ISO timestamp (`2026-10-20T03:15:42Z`). A blank or malformed value stops the run.

#### Curation status

Every output row is looked up in the cache and keeps its own `MismatchType`. `CurationStatus` is `ACKNOWLEDGED` or `CHANGED_FROM_CACHE`.

| CurationStatus | Meaning |
|----------------|---------|
| `ACKNOWLEDGED` | The latest cache row for the item has the same `MismatchType` and `Differences` |
| `CHANGED_FROM_CACHE` | The latest cache row for the item has a different `MismatchType` or `Differences` |

For later comparisons, pass the curation file as:

```bash
-curation_cache '/comparisons/pipeline_curation.tsv'
```

Example `pipeline_curation.tsv`:

```tsv
SampleId	Category	MismatchType	CurationStatus	Key	Differences	OldSource	NewSource	CurationType	Severity	Verdict	CurationComment	CuratedBy	CuratedDate
DEMO_001	GENE_COPY_NUMBER	VALUE		MYC	MinCopyNumber(6/8)	v34	v35	MANUAL	LOW	NEW_CORRECT	Expected copy number change	reviewer_a	2026-10-20
DEMO_001	GENE_COPY_NUMBER	VALUE		MYC	MinCopyNumber(8/10)	v35	v36	BATCH_INFERRED	NONE		Recorded difference	reviewer_a	2026-10-21T09:30:00Z
DEMO_002	GENE_COPY_NUMBER	VALUE		MYC	MinCopyNumber(2/5)	v35	v36	MANUAL	CRITICAL	UNCLEAR	Investigate	reviewer_b	2026-10-21
DEMO_001	RNA_FUSION	OLD_ONLY		EML4_ALK chr2:3000-chr2:4000		v35	v36	MANUAL	LOW	OLD_CORRECT	Previously reviewed filtering change	reviewer_a	2026-10-21
DEMO_001	DRIVER	NEW_ONLY		AMP_ERBB2		v35	v36	MANUAL	MEDIUM	NEW_CORRECT	Newly reported driver	reviewer_b	2026-10-21
DEMO_001	PURITY	VALUE			Purity(0.65/0.68)	v35	v36	MANUAL	LOW	BOTH_CORRECT	Small purity change	reviewer_a	2026-10-21
DEMO_001	RNA_FUSION	INVALID_NEW				v35	v36	MANUAL	NONE		Missing NEW input	reviewer_a	2026-10-21
```

For the next comparison, use these options alongside the source directories and output directory:

```bash
-sample DEMO_001 -categories GENE_COPY_NUMBER -match_level DETAILED \
-source_label_old v34 -source_label_new v37 \
-curation_cache /comparisons/pipeline_curation.tsv
```

If the v34 vs v37 comparison finds two MYC changes, `MinCopyNumber` 6 to 12 and `MaxCopyNumber` 10 to 12, both go in one row:

- `Prev*` shows the latest cache row for the item, which is the last one in the file, and `CurationStatus` compares the run's row with
  that row.
- A blank `CurationStatus` means the cache has no row for the item. Otherwise, it shows `ACKNOWLEDGED` or `CHANGED_FROM_CACHE`.
- Cache entries absent from the current comparison are not output.

The `Prev*` columns are `PrevComparison` (the cache row's `OldSource>NewSource`), `PrevMismatchType`, `PrevDifferences`,
`PrevCurationType`, `PrevSeverity`, `PrevVerdict`, `PrevCurationComment`, `PrevCuratedBy` and `PrevCuratedDate`.

| Differences | MismatchType | CurationStatus | PrevComparison | PrevDifferences |
|-------------|--------------|----------------|----------------|-----------------|
| MinCopyNumber(6/12);MaxCopyNumber(10/12) | VALUE | CHANGED_FROM_CACHE | v35>v36 | MinCopyNumber(8/10) |

A MYC row with only `MinCopyNumber(8/10)` gives `CurationStatus=ACKNOWLEDGED`. `MinCopyNumber(6/8)` gives `CHANGED_FROM_CACHE`, because
only the latest cache row counts.

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
| Filter                   | Exact                   |
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
| VariantCopyNumber        | Threshold [0.5, 15%]    |
| TumorSupportingReadCount | Threshold [1, 20%]      |
| TumorTotalReadCount      | Threshold [1, 20%]      |
| AF                       | Threshold [0.2]         |

### Germline Variant
Data key: SampleId, Chromosome, Position, Ref, Alt and VariantType (SNP/MNP/INDEL/UNDEFINED)

| Field                    | Match Type & Thresholds |
|--------------------------|-------------------------|
| Reported                 | Exact                   |
| Filter                   | Exact                   |
| Gene                     | Exact                   |
| CanonicalEffect          | Exact                   |
| CanonicalCodingEffect    | Exact                   |
| CanonicalHgvsCoding      | Exact                   |
| CanonicalHgvsProtein     | Exact                   |
| OtherReportedEffects     | Exact                   |
| Tier                     | Exact                   |
| Hotspot                  | Exact                   |
| Qual                     | Threshold [50, 20%]     |
| VariantCopyNumber        | Threshold [0.5, 15%]    |
| TumorSupportingReadCount | Threshold [1, 20%]      |
| TumorTotalReadCount      | Threshold [1, 20%]      |
| AF                       | Threshold [0.2]         |

### Germline Amp-Del
Data key: SampleId, Gene

| Field              | Match Type & Thresholds |
|--------------------|-------------------------|
| Reported           | Exact                   |
| GermlineStatus     | Exact                   |
| TumorStatus        | Exact                   |
| GermlineCopyNumber | Threshold [0.5, 15%]    |
| TumorCopyNumber    | Threshold [0.5, 15%]    |

### Drivers (Linx and Purple)
Data key: SampleId, GeneId, TranscriptId, Driver-type

| Field            | Match Type & Thresholds |
|------------------|-------------------------|
| LikelihoodMethod | Exact                   |
| Likelihood       | Threshold [0.1]         |
| MinCopyNumber    | Threshold [0.5, 15%]    |
| MaxCopyNumber    | Threshold [0.5, 15%]    |
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

| Field        | Match Type & Thresholds |
|--------------|-------------------------|
| Reported     | Exact                   |
| BreakendInfo | Exact - see below       |

A diff is raised per gene if, for any matched breakend the reported status, region type, coding type, or next splice exon rank differs between the two runs,
or a breakend/SV present in one run has no match in the other.

### Germline SV
Data key: SampleId, Gene

| Field        | Match Type & Thresholds |
|--------------|-------------------------|
| Reported     | Exact                   |
| BreakendInfo | Exact - see below       |

Compared the same way as [Disruptions](#disruptions) above.

### Cuppa
Data key: SampleId, ClassifierName, DataType

| Field           | Match Type & Thresholds |
|-----------------|-------------------------|
| TopCancerType   | Exact                   |
| Probability     | Threshold [0.1]         |

### Cuppa Image
Data key: SampleId, logical image name (`cuppa_vis_image`)

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
Data key: SampleId, genes

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
| SplicedFragmentPerc   | Threshold [0.05]        |
| UnsplicedFragmentPerc | Threshold [0.05]        |
| AltFragmentPerc       | Threshold [0.05]        |
| ChimericFragmentPerc  | Threshold [0.05]        |
| SplicedGeneCount      | Threshold [10, 1%]      |
| FragLength5th         | Threshold [5%]          |
| FragLength50th        | Threshold [5%]          |
| FragLength95th        | Threshold [5%]          |
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
| AdjTPM   | Threshold [5%]          |

### Novel Splice Junction
Data key: SampleId, GeneName, Chromosome, JunctionStart, JunctionEnd

| Field         | Match Type & Thresholds |
|---------------|-------------------------|
| Type          | Exact                   |
| FragmentCount | Threshold [5, 5%]       |
| RegionStart   | Exact                   |
| RegionEnd     | Exact                   |

### RNA Fusion
Data key: SampleId, FusionName, ChromosomeUp, PositionUp, ChromosomeDown, PositionDown

| Field        | Match Type & Thresholds |
|--------------|-------------------------|
| KnownType    | Exact                   |
| JuncTypeUp   | Exact                   |
| JuncTypeDown | Exact                   |
| SplitFrags   | Threshold [5, 5%]       |
| Filter       | Exact                   |

A fusion passing in only one source is looked up in the other source's unfiltered fusions file to show why it was filtered.

### Copy Number
Only runs in DETAILED mode.

Data key: SampleId, Chromosome, StartPosition, EndPosition

| Field                 | Match Type & Thresholds |
|-----------------------|-------------------------|
| CopyNumber            | Threshold [0.5, 15%]    |
| MajorAlleleCopyNumber | Threshold [0.5, 15%]    |
| Method                | Exact                   |

### Gene Copy Number
At REPORTABLE, only genes that Purple reports are compared.

Data key: SampleId, Gene

| Field         | Match Type & Thresholds |
|---------------|-------------------------|
| MinCopyNumber | Threshold [0.5, 15%]    |
| MaxCopyNumber | Threshold [0.5, 15%]    |
| DriverType    | Exact                   |

`MinRegionStart`, `MinRegionEnd`, `TranscriptId` and `IsCanonical` are shown in `OldValues` and `NewValues` but not compared.
