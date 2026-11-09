# Viridian

Viridian detects the presence of oncology-related viruses in samples from whole genome sequencing data.
The application functions by identifying reads with viral sequences and identifying integrations of the virus genome into the human genome (via structural variants).
Viridian runs downstream of ESVEE, Purple and BamMetrics.

Viridian was initially designed to be a replacement for the VirusBreakend [[1]](#references) and Virus Interpreter [[2]](#references) tools.

## Configuration

### Required Arguments

| Argument              | Type   | Description                                                      |
|-----------------------|--------|------------------------------------------------------------------|
| sample                | String | Sample ID.                                                       |
| tumor_bam             | Path   | Path to sample tumor BAM/CRAM, with index.                       |
| esvee_unfiltered_vcf  | Path   | Path to sample's ESVEE unfiltered VCF.                           |
| purple_dir            | Path   | Path to sample's Purple output directory.                        |
| bam_metrics_tumor_dir | Path   | Path to sample's tumor BamMetrics output directory               |
| ref_genome            | Path   | Reference genome, same as used for aligning the tumor BAM/CRAM.  |
| virus_genomes_fasta   | Path   | Path to virus genome resource FASTA.                             |
| virus_genomes_tsv     | Path   | Path to virus genome resource TSV.                               |
| virus_groups_tsv      | Path   | Path to virus group resource TSV.                                |
| output_dir            | Path   | Directory in which to store output. Created if it doesn't exist. |

### Optional Arguments

| Argument                  | Type                           | Default                     | Description                                                                             |
|---------------------------|--------------------------------|-----------------------------|-----------------------------------------------------------------------------------------|
| virus_bwa_index_image     | Path                           | `{virus_genomes_fasta}.img` | Path to BWA-MEM GATK JNI image file of the virus genome FASTA.                          |
| candidate_read_source     | String                         | `FULLY_UNMAPPED_AND_DECOY`  | Which subset of reads to search for viral reads. `FULLY_UNMAPPED_AND_DECOY` or `ALL`.   |
| extra_virus_decoy_contigs | Comma-separated list of string | (none)                      | Additional reference genome contigs which are considered virus genome.                  |
| kmer_filter               | Flag                           | (none)                      | If set, use k-mers to filter candidate viral reads before alignment (performance only). |
| bwa_lib                   | Path                           | Search in current directory | Path to BWA-MEM shared library object.                                                  |
| align_batch_size          | Integer                        | 100000                      | Number of candidate viral reads to align at once with BWA-MEM.                          |
| threads                   | Integer                        | 1                           | Number of threads to use for viral read extraction and alignment.                       |
| output_id                 | String                         | (none)                      | Additional identifier included in output file naming.                                   |
| log_level                 | String                         | `info`                      | `all`/`trace`/`debug`/`info`/`warn`/`error`/`fatal`/`off`                               |

## Example Usage

```shell
java -jar viridian.jar \
  -sample COLO829v003T \
  -tumor_bam /sample_data/COLO829v003T.bam \
  -esvee_unfiltered_vcf /sample_data/esvee/COLO829v003T.esvee.unfiltered.vcf.gz \
  -purple_dir /sample_data/purple/ \
  -bam_metrics_tumor_dir /sample_data/bam_metrics/ \
  -ref_genome /ref_data/refGenome.38.fasta \
  -virus_genomes_fasta /ref_data/viridian/virus_genomes.fasta \
  -virus_genomes_tsv /ref_data/viridian/virus_genomes.tsv \
  -virus_groups_tsv /ref_data/viridian/virus_groups.tsv \
  -threads 8 \
  -output_dir /sample_data/viridian/
```

## Methodology

### Virus Reference Data

Viridian uses a curated set of virus genomes on which its virus calling methodology is based.
The viruses are only those related to oncology, such as HPV, HBV, EBV, MCV, etc.
Viruses not curated in Viridian's reference cannot be detected by the application.

Viruses are grouped in a manner most useful for oncology, rather than focusing on the taxonomic classification, which allows simpler implementation.  
Viridian structures viruses like so:

- "Virus group". Group of viruses at the main granularity significant for reporting, e.g. HPV. Viridian calls at most one virus per group.
- "Virus". Virus subtype within a group, which is useful to identify within some groups. E.g. HPV type 16.
- "Virus genome". Specific base sequence observed in a particular virus. Each virus may have many genomes.

Viridian calls at most one virus per virus group, via selection of the virus genome most similar to what is in the tumor sample.

### Algorithm Overview

At a high level, Viridian calls viruses from two sources of information, reflecting the ways in which viruses may be present in cells:

1. Reads with virus genome sequence. Virus DNA is ambiently present in cells.
2. Virus genome integration into the host genome. The virus genome sequence (or part of) has inserted itself, or been accidentally joined into, the host genome.

Both may be significant to understanding the biology of the tumor, depending on the type of virus.

### Read Profiling

Detection of viruses from sequences of reads from the tumor BAM.
The goal is to select one virus genome per virus group - the "representative genome" - which best represents the virus sequence in the sample reads.
Reads are aligned against the representative genome to produce depth and coverage data used for reporting.

#### 1. Read extraction

Viridian scans the BAM for a subset of candidate viral reads which _may_ be attributable to a virus.

Types of fragments/reads which may be viral:

- Fully unmapped. The read and its paired mate are unaligned. Probably not human genome, could be virus genome.
- Half-mapped. One read of the pair is unmapped and its mate is aligned to the host genome. This could be a fragment spanning a virus genome integration site.
- Clipped (min 30b), if there is no supplementary alignment to a human contig. This could be a fragment spanning a virus genome integration site.
- Virus decoy contig. Some reference genomes have contigs which purposely contain virus genomes. E.g. GRCh38 with `chrEBV` contig. These are known virus sequence.

Virus decoy contigs of these genome assemblies have built-in support: GRCh38, GDC GRCh38.d1.vd1, GRCh37 hs37d5, GATK b37.
Other contigs may be specified with the `extra_virus_decoy_contigs` parameter.

In addition, these types of fragments/reads are excluded:

- Secondary and supplementary alignments.
- Duplicate fragments, as identified by REDUX.
- Unmapped by REDUX, where the original alignment was not virus decoy contig.

`candidate_read_source` parameter selects which subset of read types are considered, mostly for performance reasons:

- `ALL`: All the read types listed above. Slow on a large BAM.
- `FULLY_UNMAPPED_AND_DECOY`: Reads which are fully unmapped and reads aligned to virus decoy contigs. These types can be extracted without scanning the entire tumor BAM, which is much faster.

`FULLY_UNMAPPED_AND_DECOY` is the recommended configuration because it allows efficient read extraction which avoids reading the entire BAM.
The theoretical downside is missing the integration site fragments; however, we have not observed this to be problematic during our validation.

#### 2. Read filtering

If the `kmer_filter` argument is provided, Viridian filters the candidate viral reads extracted from the tumor BAM, before virus genome alignment, purely as a performance improvement.
The filter efficiently passes only reads which share a 19-mer exact match with any virus genome, using a bloom filter.

The performance benefit of the k-mer filter may differ between samples, depending on how many reads are identified as candidates but not actually containing viral sequence.

#### 3. Read alignment

Candidate viral reads are aligned to virus genomes using BWA-MEM [[3]](#references), with mostly default parameters and a minimum aligner score of 30.

BWA-MEM runs in `-a` mode, which produces all possible alignments for each sequence.
Multiple alignment allows intelligent selection of a single genome which best represents the sample virus overall (described in the next section).

Alignments which clip over the ends of the virus genome are discarded to reduce bias from linearisation origin selection.
Virus genomes are usually circular, but must be linearised to represent in FASTA.
The selection of the origin point may differ between genomes of the same virus group, resulting in biased coverage between genomes.
The simplest fix to make the alignment data more comparable is to drop the affected alignments.

#### 4. Representative genome selection

Within each virus group, each read "votes" for virus genomes depending on the alignment quality, and then the most voted genome is selected as the "representative genome".

Each read gives a total of one vote, which is split between virus genomes according to the edit distance of the read sequence from the genome.
As edit distance increases, the vote to that genome decreases (distributed to the other aligned genomes).
Each extra base different multiplies the genome's weight by 1/4.
Example: if a read aligned to two genomes, one perfectly and one with a mismatch, then the votes would be 0.8 and 0.2.

This voting methodology is important because genomes of the same virus group are often highly similar (upwards of 95%), such that choosing the best genome in isolation per-read may understate the virus group's overall support.
The effects of genome selection and sequencing errors would be amplified.
With voting, the genome with the highest overall similarity is more likely to be chosen.

Note that tools like VirusBreakend try to construct a consensus virus genome, which Viridian intentionally does not do.
Validation demonstrated that the Viridian read voting approach maintains high reporting concordance with VirusBreakend while reducing conceptual complexity and improving performance.

#### 5. Read assignment

Each read is assigned to the single representative genome it aligns to best, across all virus groups, so each read supports at most one virus.
The representative genome's depth and coverage are then measured from its assigned reads only, and this data flows to the virus reporting logic.

### Integration Calling

Viridian ingests insertion-like structural variants from ESVEE and aligns them to virus genomes to detect viral integrations.

Candidate variants are all those with insert sequence >= 30b, which captures novel insertions unattributed to the human genome.
In practice, most viral integrations are SGL variants, since ESVEE is not able to resolve the viral part of the sequence.

Reuse of ESVEE variants demonstrated high overall reporting concordance with VirusBreakend, while reducing integration calls revealed to be tenuous by ESVEE's methodology.

Variant insert sequences are aligned to the virus genomes with BWA-MEM.
Default parameters and a minimum aligner score of 30 are used.
If an alignment is found, the variant is considered a viral integration.

Currently, there are no other filters on the alignment quality because the possible characteristics of the integrations are largely unknown, and no filter showed value during validation.

### Reporting

As the final step, Viridian gathers each virus detected via read profiling or integrations, and marks them as "reported" or "not reported" (i.e. likely to be significant for this tumor).

To be reported, the virus group must:

- Be a reportable virus (see below), AND
  - Have at least 1 viral integration, OR
  - Be clonal, at least 0.5 copies per tumor cell

Reportable viruses are currently: HPV, HBV, EBV, HHV-8, MCV.

The clonality check uses Purple and BamMetrics data, and computes:

```text
expectedViralDepthPerCopy = tumorMeanDepth * tumorPurity / (tumorPurity * tumorPloidy + 2 * (1 - tumorPurity))
copiesPerTumorCell = virusMeanDepth / expectedViralDepthPerCopy
clonal = copiesPerTumorCell >= 0.5
```

If the Purple QC status is `FAIL_NO_TUMOR` or `FAIL_CONTAMINATION`, the clonality is not calculated and the virus can only be reported from viral integrations.
The success or failure of the clonality calculation is noted in the reporting output file `clonality_status` field.

In the reporting output, the virus displayed is the consensus of the read profiling and integration calling.
If the two paths indicate different viruses in the virus group, this is noted in the reporting output file `qcStatus` field with the value `AMBIGUOUS_VIRUS`.

## Output

Main output:

| File       | Description                                                                  |
|------------|------------------------------------------------------------------------------|
| report.tsv | Overall virus calling with reporting status. One row per called virus group. |

`report.tsv` file fields:

| Field                  | Type    | Description                                                                                     |
|------------------------|---------|-------------------------------------------------------------------------------------------------|
| name                   | String  | Human readable name of the called virus.                                                        |
| qcStatus               | String  | `NO_ABNORMALITIES` / `AMBIGUOUS_VIRUS`                                                          |
| integrations           | Integer | Count of viral integrations.                                                                    |
| interpretation         | String  | Reporting type used by hmftools downstream.                                                     |
| percentageCovered      | Number  | Percentage of representative virus genome bases aligned in reads, 0-100.                        |
| meanCoverage           | Number  | Average depth of representative virus genome.                                                   |
| reported               | Boolean | True if significant virus presence.                                                             |
| clonality_status       | String  | `CLONAL` / `NOT_CLONAL` / `NOT_PROFILED` / `PURPLE_FAIL_NO_TUMOR` / `PURPLE_FAIL_CONTAMINATION` |
| copies_per_tumor_cell  | Number  | Average number of copies of the representative virus genome in each tumor cell.                 |
| virus_group            | String  | The virus group of this call.                                                                   |

Informational/visualisation/debugging outputs:

| File                       | Description                                                                                                          |
|----------------------------|----------------------------------------------------------------------------------------------------------------------|
| candidates.fasta           | Candidate viral read sequences.                                                                                      |
| all.bam                    | Alignments of candidate viral reads to all virus genomes (BWA-MEM `-a` mode).                                        |
| representative.bam         | Alignments of candidate viral reads to their best representative virus genome only.                                  |
| genome_support.tsv         | Per-virus-genome support information from read alignment. One row per virus genome with any alignments.              |
| profile.tsv                | Result of the read profiling, mainly the representative genome info. One row per virus group with any aligned reads. |
| integration_candidates.tsv | Candidate integration structural variants with virus genome alignment info.                                          |
| manifest.tsv               | Application configuration manifest.                                                                                  |

## Known Issues

- Alignments clipping over the endpoints of the virus genome are excluded, reducing support statistics.
- Viral integrations on contigs which are not primary chromosomes are not detected since they are intentionally not called by ESVEE.
- Viral integrations with two fragments support where the second fragment soft clip is too short, will not be called, despite evidence that the fragment sequence is viral.
- Observed possibility that if multiple viral integrations occur at the same location, ESVEE may only call one.
- Observed possibility that sometimes ESVEE doesn't call a viral integration junction supported by only two fragments.
- REDUX doesn't mark duplicate fragments where the mates align to different homologous regions in `chrEBV`.
- Purple QC failure means viruses cannot be reported from clonality alone, even if the virus genome depth is high in absolute numbers.

## References

1. VirusBreakend. <https://pubmed.ncbi.nlm.nih.gov/33973999>
2. Virus Interpreter. <https://github.com/hartwigmedical/hmftools/tree/virus-interpreter-v1.7.2/virus-interpreter>
3. BWA-MEM. <https://arxiv.org/abs/1303.3997>
