# TARS

**TARS** (Transcript Alignment for RNA Splicing) makes `bwa-mem2` splice-aware for RNA reads. TARS aligns the reads against the genome with the
transcriptome added using `bwa-mem2`, then rewrites the result back to genome coordinates. The output is an ordinary genomic RNA
BAM (no transcript contigs, spliced reads carried as `N` gaps) ready for REDUX and ISOFOX.

## Contents

* [What TARS does](#what-tars-does)
* [How to Run TARS](#how-to-run-tars)
* [What a read goes through](#what-a-read-goes-through)
    * [Step 0: Translate transcriptome alignments to reference genome](#step-0-translate-transcriptome-alignments-to-reference-genome)
    * [Step 1: Score short overhangs against the reference genome, collapse weak scoring ones](#step-1-score-short-overhangs-against-the-reference-genome-collapse-weak-scoring-ones)
    * [Step 2: Resolve supplementary records into splice junction candidates](#step-2-resolve-supplementary-records-into-splice-junction-candidates)
        * [Step 2.1: Pick the main alignment for a supplementary record](#step-21-pick-the-main-alignment-for-a-supplementary-record)
        * [Step 2.2a: Resolve supplementary records into splice junctions](#step-22a-resolve-supplementary-records-into-splice-junctions)
        * [Step 2.2b: On a successful resolve](#step-22b-on-a-successful-resolve)
    * [Step 3: Decide which alignments to keep for a read](#step-3-decide-which-alignments-to-keep-for-a-read)
    * [Step 4: Emit records](#step-4-emit-records)

## What TARS does

A normal genome aligner like `bwa-mem2` matches a read against one continuous stretch of genome. It has no idea where
introns are, so a read that jumps over one gets cut short or placed in the wrong spot.

TARS fixes this without changing the aligner:

1. **`SpliceFastaBuilder`** (run once per Ensembl release) concatenates each multi-exon transcript's exon sequences into
   a transcript contig (`*_tx`, introns removed); these contigs are the transcriptome. It also writes a sidecar TSV
   mapping each contig's intervals back to their genomic exon spans.
2. Append the transcriptome to the genome FASTA and index with `bwa-mem2`.
3. Align the RNA reads with `bwa-mem2` as usual. A read that jumps an intron now has a continuous place to land.
4. **`TarsApplication`** lifts each read back to its real genome position, marks the skipped intron as a gap (`N`), fixes
   things up (tags, mate info, confidence), and writes a sorted, indexed BAM.
5. Feed the new splice-aware records BAM file into REDUX (dedup), then ISOFOX.

![TARS pipeline](doc/tars.svg)

## How to Run TARS

### Build the transcript reference (SpliceFastaBuilder)

```
java -cp tars.jar com.hartwig.hmftools.tars.fasta.SpliceFastaBuilder
    -ensembl_data_dir /ref_data/ensembl_data_cache/38/
    -ref_genome /path_to_fasta/genome.fasta
    -ref_genome_version V38
    -output_dir /path_to_output/
```

Two files are written:

* `ref_genome_v38_rna_contigs.fasta` - the transcript contigs.
* `ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv` - the contig sidecar (intervals mapped back to genomic exon spans).

Concatenate the FASTA onto the genome FASTA and `bwa-mem2 index` the result before aligning.

### Run TARS (TarsApplication)

```
java -jar tars.jar
    -sample COLO829T
    -input_bam COLO829T.lane_01.bwa_tx.namegrouped.bam,COLO829T.lane_02.bwa_tx.namegrouped.bam
    -ref_genome /path_to_fasta/genome_plus_tx.fasta
    -contig_sidecar /path_to/ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv
    -rna_unmap_regions /ref_data/rna/38/rna_excluded_regions.38.tsv
    -bamtool /path/to/samtools
    -output_dir /path_to_output/
    -threads 24
```

### Output files

Every file is named `<sample>.tars.<...>`. Two are written by default:

* `<sample>.tars.bam` (+ `.bai`) - the lifted, coordinate-sorted genomic BAM, ready for REDUX.
* `<sample>.tars.summary.tsv` - a counts summary of what liftback did.

`-output_id chr1_slice` inserts the token into every name: `<sample>.tars.chr1_slice.bam`.

### Flags

**Required**

| Flag               | Description                                                                  |
|--------------------|------------------------------------------------------------------------------|
| sample             | Sample ID; prefix for each output file (`<sample>.tars.*`)                  |
| input_bam          | bwa-mem2 output against the combined FASTA, **name-grouped** (not coord-sorted). Separate with `,` for multiple lane BAMs |
| ref_genome         | The same combined genome + transcript FASTA used at alignment                 |
| contig_sidecar     | Contig sidecar TSV from `SpliceFastaBuilder` (`*.rna_contigs_mappings.tsv`)    |
| bamtool            | Path to samtools or sambamba; concatenates, sorts, and indexes the output      |
| output_dir         | Directory for the lifted BAM and summary file                                     |

**Optional**

| Flag               | Default | Description                                                              |
|--------------------|---------|--------------------------------------------------------------------------|
| output_id          | (none)  | ID inserted into every output filename |
| rna_unmap_regions  | (none)  | Curated excluded regions (rRNA / 7SL / multi-map zones) whose reads are unmapped in the lifted output using REDUX SAM conventions; see [rna_excluded_regions.38.tsv](https://source.cloud.google.com/hmf-pipeline-development/common-resources-public/+/master:rna/38/rna_excluded_regions.38.tsv) |
| threads            | 1       | Worker threads; reads process in parallel per read-group |

**Tuning thresholds**

| Flag                            | Default   | Description                                                   |
|---------------------------------|-----------|---------------------------------------------------------------|
| supp_implied_min_intron_length  | 21        | Min implied intron length for a primary+supp merge            |
| supp_implied_max_intron_length  | 1000000   | Max implied intron length for a primary+supp merge            |

Note: no `ensembl_data_dir` - liftback reads exon/junction annotation from the sidecar (only `SpliceFastaBuilder` needs
Ensembl).

### Upstream bwa-mem2 flags

Not tars config, but liftback depends on them.

| Setting | Value | What it is |
|---|---|---|
| `-T` | 19 | bwa-mem2 minimum alignment score to output; set below the default 30 to retain short-anchor supplementaries for Step 2 |
| `-h` | 75 | bwa-mem2 XA hit cap used by the TARS alignment |

## What a read goes through

After bwa-mem2 alignment, a read that spans an exon boundary or has supplementaries around novel junctions is processed by
TARS through these steps in order:

```
Step 0  Translate   lift every primary and supplementary placement, including XA alts
Step 1  Overhang    re-evaluate each overhang and collapse the weak scoring ones
Step 2  Merge       pick one alignment per supplementary record, then try splice merges
Step 3  Decide      choose the mate placements together, then keep the remaining XA alternates
Step 4  Emit        update the records and write the BAM
```

### Step 0: Translate transcriptome alignments to reference genome

Every read's transcriptome alignment is translated to genomic coordinates, with introns re-inserted as `N` gaps / splice junctions.

![translate the read to the genome](doc/translate.svg)

### Step 1: Score short overhangs against the reference genome, collapse weak scoring ones

A short overhang (`<= 12M`) next to a splice junction at a read end is re-scored using bwa-mem2-style scoring against the
reference genome. There are 3 cases:

**1a.** With a soft clip: keep the junction if the overhang scores > 5; otherwise drop the `N` junction and walk the soft
clip onto the reference genome, leaving a contiguous alignment.

![1 splice junction](doc/overhang_one_junction.svg)

**1b.** With multiple splice junctions: keep the junction if the short overhang aligns positively (AS > 0), otherwise collapse
it only when the intronic reference AS > short overhang AS.

![more than 1 splice junction](doc/overhang_two_junctions.svg)

**1c.** With no soft clip and a single junction: not checked, no intervention.

### Step 2: Resolve supplementary records into splice junction candidates

`bwa-mem2` is run with `-T 19`, allowing short-anchor supplementary alignments at junction sites (annotated or novel) to
be kept.

#### Step 2.1: Pick the main alignment for a supplementary record

For `MAPQ > 0`, TARS picks its sole alignment. For `MAPQ 0`, TARS picks one alignment (main + `XA` entries) from each
supplementary record (within 1 Mb):

1. shortest `DEL`
2. shortest `DUP`
3. shortest `INV`
4. deterministic random

A selected `XA` alignment becomes the supplementary's main alignment, whether or not a merge succeeds.

#### Step 2.2a: Resolve supplementary records into splice junctions

TARS tries the selected alignments against the primary splice chain. A merge requires:

- the same chromosome and strand
- simple CIGARs with complementary terminal soft clips
- full read coverage, at most 5 bp overlap, and no reference overlap
- the implied intron length is within [`supp_implied_min_intron_length`, `supp_implied_max_intron_length`]
- at most one supplementary merge per terminal soft clip

A successful merge becomes a candidate for Step 3. If it wins, the absorbed supplementary records are dropped.

![merge supplementary record to primary](doc/rescue_via_supplementary.svg)

#### Step 2.2b: On a successful resolve

The junction position may still be ambiguous. TARS uses this order:

1. an annotated junction (Ensembl)
2. a canonical `GT-AG` splice motif, then semi-canonical
3. the mate's already-resolved junction
4. the midpoint of the ambiguous read range, rounded down

Ties at the chosen annotated or motif tier use a deterministic read-seeded choice.

### Step 3: Decide which alignments to keep for a read

A read now has its own alignment plus any `XA` alternate alignments: each a genomic (ref) or translated transcriptome
(tx) alignment, plus any supplementary-supported merge candidates. TARS picks one as the primary, keeps only the relevant
`XA`, and sets its `MAPQ`. Every read lands in one of three buckets:

- **B1. Ref only:** keep BWA's primary. A genomic `MAPQ 0` read with no `XA` is treated as over-cap and unmapped.

- **B2. Ref/tx agreement:** ref and tx alignments lift to the same contiguous locus and CIGAR; keep BWA's primary.

- **B3. Multi-mapper:** for a `MAPQ 0` pair, TARS selects both mates together (within 1 Mb; distances exclude introns):
    1. closest mates under 1 kb apart (`9M2209N142M` measures from its `142M` exon, not across the intron)
    2. shortest `DEL`
    3. shortest `DUP`
    4. shortest `INV`
    5. highest combined alignment score
    6. deterministic random

  A mate with positive MAPQ is not moved; it anchors the pair. Supplementary merges from Step 2 still apply at any MAPQ.

### Step 4: Emit records

TARS sets the MAPQ in this order:

1. keep a positive input MAPQ
2. raise `MAPQ 0` to 60 at a single locus, unless the pick was a random tie
3. leave an `XS == AS` tie at 0 unless the winner came from a transcript contig or an annotated exon
4. take `max(primary, supplementary)` for a supplementary merge, or 60 at a single locus

TARS then finalises each record:

- write the remaining alignments to `XA`, rebuild `SA`, update `AS`, `XS`, `NH`, `NM` and mate fields, remove `MD`
- drop absorbed, duplicate, unliftable, excluded, `AS < 30`, and supplementaries with no surviving `SA` partner
- unmap unliftable, excluded, over-cap, and `AS < 30` primaries using REDUX conventions
- skip expected `inter-transcript spacer` misses, log other lift failures
