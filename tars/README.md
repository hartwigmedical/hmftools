# TARS

**TARS** (Transcript Alignment for RNA Splicing) makes `bwa-mem2` splice-aware for RNA reads. TARS aligns the reads against the genome with the
transcriptome added using `bwa-mem2`, then rewrites the result back to genome coordinates. The output is an ordinary genomic RNA
BAM (no transcript contigs, spliced reads carried as `N` gaps) ready for REDUX and ISOFOX.

## Contents

* [What TARS does](#what-tars-does)
* [How to Run TARS](#how-to-run-tars)
* [What a read goes through](#what-a-read-goes-through)
* [What TARS changes on a read](#what-tars-changes-on-a-read)
* [In-depth technical implementation](#in-depth-technical-implementation)

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

![TARS pipeline](doc/tars.png)

## How to Run TARS

<details>
<summary><strong>Commands, output files and flags</strong></summary>

### Build the transcript reference (SpliceFastaBuilder)

```
java -cp tars.jar com.hartwig.hmftools.tars.fasta.SpliceFastaBuilder
    -ensembl_data_dir /ref_data/ensembl_data_cache/38/
    -ref_genome /path_to_fasta/genome.fasta
    -ref_genome_version V38
    -output_dir /path_to_output/
```

Outputs `ref_genome_v38_rna_contigs.fasta` and `ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv`. Concatenate the
FASTA onto the genome FASTA and `bwa-mem2 index` the result before aligning.

### Run TARS (TarsApplication)

```
java -jar tars.jar
    -sample ACTN01020030T
    -input_bam ACTN01020030T.bwa_tx.namegrouped.bam
    -ref_genome /path_to_fasta/genome_plus_tx.fasta
    -ref_genome_version V38
    -contig_sidecar /path_to/ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv
    -rna_unmap_regions /ref_data/rna/38/rna_excluded_regions.38.tsv
    -bamtool /path_to_samtools/
    -output_dir /path_to_output/
    -threads 24
```

### Output files

Named `<sample>.tars[.<output_id>][.<stage>].ext`. You get two things:

* `ACTN01020030T.tars.bam` (+ `.bai`) - the lifted, coord-sorted genomic BAM, ready for REDUX.
* `ACTN01020030T.tars.summary.tsv` - a counts summary of what liftback did.

With `-output_id chr1_slice` the token is inserted into every name: `ACTN01020030T.tars.chr1_slice.bam`. The per-record
debug TSVs are written only with `-write_liftback_tsv` (a ~100GB file; only for per-read detail the summary can't give).

### Flags

**Required**

| Flag               | Description                                                                  |
|--------------------|------------------------------------------------------------------------------|
| sample             | Sample ID. prefix to each output file (`<sample>.tars.*`)                  |
| input_bam          | bwa-mem2 output against the combined FASTA, **name-grouped** (not coord-sorted)|
| ref_genome         | The same combined genome + transcript FASTA used at alignment                 |
| ref_genome_version | `V37` or `V38`                                                                |
| contig_sidecar     | Contig sidecar TSV from `SpliceFastaBuilder` (`*.rna_contigs_mappings.tsv`)    |
| bamtool            | samtools path (used to decompress the input and sort + index output)          |
| output_dir         | Directory for the lifted BAM and summary                                      |

**Optional**

| Flag               | Default | Description                                                              |
|--------------------|---------|--------------------------------------------------------------------------|
| output_id          | (none)  | id inserted into every output name, e.g. `-output_id chr1_slice` -> `ACTN01020030T.tars.chr1_slice.bam` |
| rna_unmap_regions  | (none)  | TSV of curated regions (rRNA / 7SL / multi-map); reads lifting in are unmapped. TODO: link file |
| write_liftback_tsv | off     | Per-record debug TSVs; off by default (creates a ~100GB file)            |
| threads            | 1       | Worker threads; reads process in parallel per read-group (prod uses 32 threads) |
| log_level          | INFO    | `-log_debug` for DEBUG, or `-log_level DEBUG_2` for per-read liftback detail |

**Tuning thresholds**

| Flag                        | Default   | Description                                                       |
|-----------------------------|-----------|------------------------------------------------------------------|
| rescue_min_intron           | 21        | Min intron length for a primary+supp merge                       |
| rescue_max_intron           | 1000000   | Max intron length for a primary+supp merge                       |
| rescue_min_anchor_overhang  | 3         | Min matched bases each side of a merged junction                 |
| rescue_max_chain_depth      | 2         | Max supplementaries merged into one read (one per side)          |
| rescue_softclip_tolerance   | 5         | Primary/supp overlap tolerance when snapping to a junction       |
| rescue_max_boundary_shift   | 8         | Max over-extended bases trimmed when probing a boundary          |
| rescue_min_partial_match_run| 11        | Min exon-proximal matched run for a partial ref-verify rescue    |
| tail_min_extension          | 3         | Min softclip length considered, and min ref-matching bases converted to M |
| tail_max_extension          | 30        | Max bases walked into a softclip                                 |

Note: no `ensembl_data_dir` - liftback reads exon/junction annotation from the sidecar (only `SpliceFastaBuilder` needs
ensembl). What each threshold trades off if you change it is in
[tuning thresholds and their effect](#in-depth-technical-implementation) below.

</details>

## What a read goes through

Liftback translates each tx-contig alignment to genomic coordinates, reconciles every candidate cigar so the ref-vs-tx
choice is made, picks the primary, then refines that primary. Execution order (each pass feeds the next):

```
translate every candidate (self + lifted XA alts)                         [Step 0]
  -> reconcile each candidate: collapse, then tail-extend                 [Steps 1 + 2]
  -> discriminate ref vs tx on the reconciled cigars                      [Step 3]
  -> rescue via supplementary: merge the primary with a supp              [Step 4]
  -> canonicalize: slide the junction onto a splice motif                 [Step 5]
```

### Step 0: translate to the genome

`ContigTranslator` rewrites al transcript-contig alignments (including XA alts) to the reference genome for comparision 

- A read that crosses an exon boundary becomes splice aware
- A read contained within one exon crosses no boundary and lifts unchanged

![translate the read to the genome](doc/translate.svg)

### Step 1: re-check short junctions

A junction is only trusted with at least `MIN_JUNCTION_ANCHOR` (8 bp) matched. Any junction with a shorter anchor is
re-scored against the reference genome. The higher score wins:

- if the contiguous extension scores at least as well, the junction is dropped and the matching bases become `M`
- only the best-scoring prefix is reclaimed, so a mismatched tail stays softclipped (e.g. last 2 of 5 mismatch give `149M 2S`, not `151M`)
- if the split scores strictly better, the junction is kept

![re-check short junctions](doc/terminal_micro_junction_collapse.svg)

### Step 2: re-check softclips within threshold

A softclip of at least `tail_min_extension` is walked into the contiguous genome (no gap) and re-scored.
The best-scoring prefix, up to `tail_max_extension`, is turned from softclip (`S`) into matched (`M`):

- the best prefix must be at least `tail_min_extension`, otherwise the tail is left clipped
- only the matching prefix is converted; any mismatched remainder stays softclipped (e.g. `131M 20S` becomes `139M 12S`)
- if the walk stalls right where an annotated intron begins, the tail is left for rescue (Step 4) instead

![re-check softclips within threshold](doc/tail_extension.svg)

### Step 3: decide which alignment to keep

A read can match the genome, the transcriptome, or both. When more than one alignment is present, TARS makes a decision
on which one to keep. On a confident placement it updates the record's CIGAR, pos & tags; otherwise the read is left untouched.

The three most common decisions:

**Transcript only, the core splice**

The read matched only the transcriptome, so bwa scored it MAPQ 0 (it could not tell the copies apart). The copies all
point to one genomic spot, so TARS places it there with the intron as an `N` gap, drops the duplicates, and marks it confident.

![transcript only, the core splice](doc/core_splice_tx_only.svg)

**Both copies agree**

The read matched the genome and the transcriptome at the same place with the same shape (no intron). TARS keeps the
genome placement and drops the duplicate transcriptome alignment.

![both copies agree](doc/concordant.svg)

**Transcript crosses a junction, genome lies flat**

The transcriptome alignment carries an `N` junction; a genomic alignment at another locus is contiguous (no `N`). bwa
often makes the contiguous one the primary to keep the read near its mate (proper-pair rescue), even though the spliced
placement is correct. TARS keeps the junction-spanning placement: it drops the contiguous alt, or swaps to the spliced
placement when bwa had made the contiguous one primary.

![transcript crosses a junction, genome lies flat](doc/swap_junction_over_contiguous.svg)

The remaining cases:

- **Genome only** - matched only the genome; nothing to decide, placement untouched.
- **Several spots, no winner** - matches multiple places with no clean splice signal; left multi-mapped (MAPQ 0).
- **No rule fits** - one spot with both copies but a shape no rule covers; bwa's choice is kept.
- **Genome reads through a transcript gap (rare)** - the transcript opens a tiny gap, but the genome reads straight
  through; the genome placement wins.

### Step 4: merge supplementary record to primary

bwa-mem2 is run with a low alignment-score floor (`-T 19`) to keep short-anchor supplementary alignments at junction
sites (annotated or novel). TARS merges such a supplementary back into its primary across the intron, giving one spliced 
primary with an `N` gap, and drops the supplementary. The merge requires:

- the primary and supplementary cover complementary halves of the read
- the implied intron is within [`rescue_min_intron`, `rescue_max_intron`]
- at least `rescue_min_anchor_overhang` matched each side of the junction
- the junction snaps to a nearby annotated junction or splice motif when one exists; otherwise it merges at bwa's split
  point, so novel junctions are merged too

The merged primary's MAPQ is capped at `SUPP_RESCUE_MAPQ_CAP` (55), marking a constructed alignment.

![merge supplementary record to primary](doc/rescue_via_supplementary.svg)

### Step 5: slide the junction onto a splice motif

A junction can sit a base or two off the splice motif. TARS slides the intron up to 5 bp onto the best motif, taking the
smallest shift where every moved base still matches.

![slide the junction onto a splice motif](doc/canonicalize.svg)

## What TARS changes on a read

For each read, TARS rewrites these fields.

**RNAME / POS / strand / CIGAR**

- set to the resolved genomic placement, the output of the six steps above, with an `N` gap per intron
- cleared to `*` / `0` / `*` when the read is unmapped

**MAPQ**

- raised to 60 when the read clearly belongs at one spot: the discriminator picked a transcript placement; or a single
  best-scoring locus bwa scored 0 (`numLoci == 1`, `inputMapq == 0`) with no hidden tie; or a transcript supplementary
  bwa scored 0
- capped at 55 after a primary+supp merge (marks a stitched alignment)
- left at bwa's value when there is no transcript evidence, or on a hidden tie (`XS == AS`, no transcript, outside any exon)
- set to 0 when the read is unmapped
- never touched by the cleanup passes (collapse, tail-extend, canonicalize)

**XA**

- rebuilt from scratch, format `chrom,±pos,CIGAR,NM;`
- an alt is kept only if it is not the chosen placement, not dropped, and does not overlap the read's own span
- cleared when the read is unmapped

**SA**

- each entry lifted to genome coords; its NM and MAPQ carried over unchanged
- an entry dropped if it fails to lift, belongs to a dropped split piece, or duplicates one already kept
- a split piece whose entries all fail to lift is dropped

**MC and mate fields**

- `MC` = the mate's lifted CIGAR
- mate RNAME / POS / strand = the mate's lifted placement
- `TLEN` = the signed 5'-to-5' distance
- an unmapped mate is parked on this read's own coordinates

**NM / MD**

- `NM` recomputed against the genome over the `M` blocks only
- `MD` always dropped, never rebuilt

**NH**

- `NH = max(numLoci, 1)`, counting distinct genomic loci among the best-scoring alignments only
- one junction shared by many transcripts counts once

**XS:A**

- set to `+` / `-` only when the lifted CIGAR has an `N` and the transcript strand is known
- bwa's `XS:i:` (a different value on the same tag name) is cleared first
- Isofox reads it for stranded junctions

**Unmapped or dropped** when the read cannot be trusted at any one place:

- its lift failed (a position outside every transcript)
- it is past bwa's locus cap (`-h 75`: MAPQ 0 with no XA)
- it scores below the floor (`AS < 30`)
- it lifts into a curated rRNA / contamination region
- a split piece whose own lift failed is moved onto its primary's coordinates instead of being unmapped

## In-depth technical implementation

Everything below is the detailed mechanics.

<details>
<summary><strong>Choosing ref vs transcript: the discriminator</strong></summary>

bwa-mem2 aligns each read to a ref contig, a transcript (`tx`) contig, or both, and picks one as primary. The
discriminator decides what to do with that primary. Every read lands in one of four **record actions** (by frequency):

- **Kept as-is** - bwa's primary is already right; nothing moves.
- **Losing alt dropped** - a ref and a tx alignment disagree; the loser is deleted from the read's alts.
- **Primary swapped** - bwa had made the loser its primary; the winner is moved into the record's main fields.
- **Left multimapped** - no rule resolves it; bwa's primary and its alts are kept untouched.

It runs on **reconciled** candidates: every alignment (the bwa primary plus each lifted XA alt) is first put through the
collapse and tail-extend passes, leaving one of three tested shapes - a contiguous full match, a trusted `N` junction,
or an unreclaimable residual softclip. Its `RefSoftClipped` / `RefFullMatch` / `TxHasNCigar` inputs are therefore
measured facts, not raw-bwa guesses: a ref softclip the genome explains contiguously is reclaimed first, so it no longer
masquerades as evidence for a tx junction. An `N` counts as a real splice only with >= 8 matched (`M`) bases each side
(`MIN_JUNCTION_ANCHOR`) - one strong junction trusts the whole read.

The plain-language walkthrough of the common cases is in [Step 3](#step-3-decide-which-alignment-to-keep); the formal
rule for every case is the decision tree below, whose `DecidingFeature` leaves (`SOLE_REF`, `JUNCTION`, ...) name each case.


#### The decision tree

Each leaf is an `(Outcome, DecidingFeature)` pair; the record action (kept / dropped / swapped / multimapped) follows
from the Outcome.
```
loci >= 2 ?
|
+-- YES (multi-locus)
|     hasRef AND hasTx ?
|     |  YES:  TxHasNCigar AND NOT RefHasNCigar -> (TX,         JUNCTION_OVER_CONTIGUOUS)  [tx faithful; ref alts contiguous]
|     |        else                             -> (UNRESOLVED, MULTIMAPPER)             [genuine multimapper, no swap]
|     |  tx only                                -> (TX,         SOLE_TX)
|     |  ref only                               -> (REF,        SOLE_REF)
|
+-- NO (single locus)
      only ref                                  -> (REF,        SOLE_REF)
      only tx                                   -> (TX,         SOLE_TX)
      both ref and tx:
         same CIGAR and no N                    -> (REF,        CONCORDANT)         [no contest]
         TxHasNCigar AND RefSoftClipped         -> (TX,         JUNCTION)           *** TX WINS ***
         TxHasNCigar AND RefFullMatch           -> (REF,        REF_READS_THROUGH)  *** REF WINS ***
         TxSoftClipAtBoundary AND RefFullMatch  -> (REF,        REF_READS_THROUGH)  *** REF WINS ***
         otherwise                              -> (UNRESOLVED, AMBIGUOUS)          [no rule fires, no swap]

  TxHasNCigar = >=1 N with >= 8 matched bases each side    RefFullMatch = ref alt with no softclip
```
Both single-locus REF-wins scenarios - the tx alignment carrying an `N` junction, and the tx alignment only softclipping
at the boundary - resolve to `(REF, REF_READS_THROUGH)`; the `TxHasNCigar` / `TxSoftClipAtBoundary` flags distinguish them.

#### Swap vs drop

Swap is rare - every contested category is drop-dominated (whole sample; drop = records - swaps):

* `JUNCTION`: 2,779,734 drop / 20,456 swap
* `JUNCTION_OVER_CONTIGUOUS`: 7,393,896 drop / 434,900 swap
* `REF_READS_THROUGH`: 17,572 drop / 12,258 swap

`REF_READS_THROUGH` (the rarest) swaps nearly as often as it drops; the others almost always just drop.

</details>

<details>
<summary><strong>Worked BAM records (before / after)</strong></summary>

Real records from a production run (sample ACTN01020030T), shown before liftback (bwa-mem2 against the combined
genome + transcript-contig FASTA) and after. Each `raw` block is the SAM record trimmed to the columns that carry the
story; QNAME, RNEXT/PNEXT/TLEN, SEQ/QUAL, MD and RG are dropped and long `XA` lists truncated. A `*_tx` reference name
means the read aligned to a transcript contig; a plain `chrN` is genomic. Column order in every block:

    FLAG  RNAME  POS  MAPQ  CIGAR  tags (NM / AS / XS / NH / XA / SA)

**SOLE_TX - core splice, MAPQ 0 -> 60**

```
chr22_tx:8725640  151M  MAPQ 0      ->      chr22:45923222  114M7762N37M  MAPQ 60
```

The read aligns flat against one of three packed transcript contigs, all tied at MAPQ 0. Liftback translates the match
across the 7,762 bp intron, the three copies collapse to a single genomic locus, the redundant alts are dropped, and
MAPQ is rescued from 0 to 60.

<details><summary>raw</summary>

```
        FLAG  RNAME     POS       MAPQ  CIGAR         tags
before  83    chr22_tx  8725640   0     151M          NM:1 AS:146 XS:146  XA:chr22_tx,-8731128,151M,1; chr22_tx,-8728343,151M,1
after   83    chr22     45923222  60    114M7762N37M  NM:1 AS:146 NH:1    XS:A:-
```
</details>

**Swap (JUNCTION_OVER_CONTIGUOUS) - MAPQ 0 -> 60**

```
chrX:54809862  151M (contiguous)  MAPQ 0      ->      chrX:54809808  10M54N141M  MAPQ 60
```

bwa places the read as a flat, contiguous 151M (MAPQ 0, tied across transcript contigs), blind to the splice. Liftback
recognises the 54 bp junction, swaps to the spliced placement, drops the redundant alts, and rescues MAPQ.

<details><summary>raw</summary>

```
        FLAG  RNAME  POS        MAPQ  CIGAR        tags
before  147   chrX   54809862   0     151M         NM:1 AS:146 XS:146  XA:chrX_tx,-7027849,151M,1; +7 more
after   147   chrX   54809808   60    10M54N141M   NM:2 AS:146 NH:1    XS:A:+
```
</details>

**Ref wins - tx alt dropped, MAPQ unchanged (54)**

```
chr1:155188564  151M  MAPQ 54   + tx alt chr1_tx 109M42S      ->      chr1:155188564  151M  MAPQ 54   (tx alt dropped; placement unchanged)
```

The genome reads through contiguously (3 mismatches); the transcript alt only anchors 109M before softclipping, never
spanning a real junction. Ref is kept, the tx alt dropped, and MAPQ is carried through unchanged - this was never a
MAPQ-0 tie, so no rescue applies.

<details><summary>raw</summary>

```
        FLAG  RNAME  POS         MAPQ  CIGAR  tags
before  163   chr1   155188564   54    151M   NM:3 AS:136 XS:109  XA:chr1_tx,+28972326,109M42S,0
after   163   chr1   155188564   54    151M   NM:3 AS:136 NH:1
```
</details>

**Rescue via supplementary - MAPQ 60 -> 55 (capped)**

```
chr22:45885973 32S113M6S (primary)  +  chr22:45885775 35M116S (supp)      ->      chr22:45885775  35M166N110M6S  MAPQ 55
```

bwa emits the read as a primary plus a short-anchor supplementary. Rescue merges them across the 166 bp intron into one
spliced primary and drops the supplementary. MAPQ is capped at 55 to flag a constructed alignment, not a bwa-native one.

<details><summary>raw</summary>

```
        FLAG  RNAME  POS        MAPQ  CIGAR          tags
before  97    chr22  45885973   60    32S113M6S      NM:0 AS:113  SA:chr22,45885775,+,35M116S,60,0       (primary)
before  2145  chr22  45885775   60    35M116S        NM:0 AS:35   SA:chr22,45885973,+,32S113M6S,60,0     (supp, FLAG 0x800)
after   97    chr22  45885775   55    35M166N110M6S  NM:0 AS:113 NH:1                                    (supp merged in, then dropped)
```
</details>

**Terminal micro-junction collapse - CIGAR only, MAPQ set by the discriminator**

```
chr22_tx:8732934  151M   --translate-->   chr22:45931220  150M18550N1M   --collapse-->   151M   MAPQ 60
```

Translation places a 1 bp anchor across an 18,550 bp "junction". The contiguous genome explains the tail at least as
well, so collapse reclaims it to a plain 151M. (MAPQ here moved 0 -> 60 via the single-locus rescue; the collapse pass
itself does not touch MAPQ.)

<details><summary>raw</summary>

```
        FLAG  RNAME     POS        MAPQ  CIGAR  tags
before  163   chr22_tx  8732934    0     151M   NM:0 AS:151 XS:151  XA: ...; chr22,+45931220,151M,0; ...
        intermediate: translated to 150M18550N1M, collapse reclaims it -> 151M   (DEBUG_2)
after   163   chr22     45931220   60    151M   NM:0 AS:151 NH:1
```
</details>

**Tail extension - CIGAR only, MAPQ unchanged**

```
chr16:22503853  27S116M8S  MAPQ 0      ->      chr16:22503826  143M8S  MAPQ 0
```

The 27 bp leading softclip matches contiguous genome, so it is walked back into the alignment (27S + 116M -> 143M);
the 8 bp trailing residual stays clipped. The read is a multimapper, so MAPQ stays 0 - tail extension never changes it.

<details><summary>raw</summary>

```
        FLAG  RNAME  POS        MAPQ  CIGAR      tags
before  99    chr16  22503853   0     27S116M8S  NM:13 AS:51 XS:53  XA: +38 entries
after   99    chr16  22503826   0     143M8S     NM:14 AS:51 NH:1   XA: +38 entries
```
</details>

**Canonicalize - CIGAR only, MAPQ set by the discriminator**

```
chr22_tx:8728628  100M4D51M   --translate-->   100M18554N51M   --canon-->   chr22:45931269  101M18554N50M  MAPQ 60
```

The boundary `4D` (a transcript-FASTA artefact) is absorbed into the intron during translation, then the junction
slides 1 bp onto the canonical splice motif. `NM` recomputes to 0 against the genome (the deletion was not real).

<details><summary>raw</summary>

```
        FLAG  RNAME     POS        MAPQ  CIGAR          tags
before  83    chr22_tx  8728628    0     100M4D51M      NM:4 AS:141 XS:141  XA:chr22_tx,-8725925,100M4D51M,4; ...
        intermediate: translated to 100M18554N51M, canonicalize slides -> 101M18554N50M   (DEBUG_2)
after   83    chr22     45931269   60    101M18554N50M  NM:0 AS:141 NH:1   XS:A:-
```
</details>

**MAPQ summary**

| scenario                          | MAPQ before | MAPQ after | why                                         |
|-----------------------------------|-------------|------------|---------------------------------------------|
| SOLE_TX                           | 0           | 60         | transcript contigs collapse to one genomic locus -> rescued |
| Swap (JUNCTION_OVER_CONTIGUOUS)     | 0           | 60         | the transcript junction placement is chosen -> rescued |
| Ref wins                          | 54          | 54         | ref was already primary, not a MAPQ-0 tie -> carried through |
| Rescue via supplementary          | 60          | 55         | a constructed primary+supp merge, capped    |
| Tail / Collapse / Canonicalize    | unchanged   | unchanged  | CIGAR-only passes; MAPQ is the discriminator's |

In short: 60 means liftback is asserting a rescue, 55 means a constructed merge, 0 means unmapped or genuinely
ambiguous, and any other value is bwa's own MAPQ carried through.

</details>

<details>
<summary><strong>Thresholds and constants</strong></summary>

All live in `TarsConstants`; the rescue and tail values are also tunable via the [flags](#flags) above. Each entry gives
the default, what it gates, and the effect of changing it. Raising a floor is the conservative direction; lowering it is
aggressive (more changes, more coincidental-match risk).

Junction rescue (primary + supp merge):

- `rescue_min_intron` / `DEFAULT_MIN_INTRON_LENGTH` (21): shortest intron a merge may bridge. Higher rejects more short
  gaps as introns (won't call a deletion a splice); lower merges small gaps, risking false junctions.
- `rescue_max_intron` / `DEFAULT_MAX_INTRON_LENGTH` (1,000,000): longest intron a merge may bridge. Higher allows merges
  across larger gaps (risk: joining unrelated loci); lower drops genuine long-intron splices.
- `rescue_min_anchor_overhang` / `DEFAULT_MIN_ANCHOR_OVERHANG` (3): matched bases kept each side of a merged junction.
  Higher = fewer chance merges but loses short-overhang splices; lower lets 1-2 bp anchors merge (~1 in 16 by chance).
- `rescue_max_chain_depth` / `DEFAULT_MAX_CHAIN_DEPTH` (2): max supps merged into one read (one per side). Higher chains
  3+ supps (recovers many-exon reads, more mis-merge risk); 1 leaves a fully-clipped middle exon unmerged.
- `rescue_softclip_tolerance` (5): primary/supp overlap tolerated when snapping to a junction. Higher tolerates bigger
  overlap (risk: double-covered bases); lower rejects the common 1-5 bp boundary disagreement and cuts yield.
- `rescue_max_boundary_shift` (8): over-extended bases trimmed when probing a boundary. Higher searches further back for
  an annotated junction (risk: a wrong nearby boundary); lower probes only bwa's exact boundary.
- `rescue_min_partial_match_run` (11): min exon-proximal matched run for a partial ref-verify rescue. Higher is safer;
  lower rescues short anchors with long unexplained tails (more false junctions).

Tail extension:

- `tail_min_extension` / `DEFAULT_MIN_EXTENSION` (3): min softclip length considered and min bases reclaimed. Higher
  ignores short clips and short matches; lower reclaims 1-2 bp (often coincidental).
- `tail_max_extension` / `DEFAULT_MAX_EXTENSION` (30): max bases walked into a softclip. Higher reclaims longer tails
  (risk: walking across a real unannotated splice); lower leaves long matching tails partly clipped.

Junctions and splice motifs (constants):

- `MIN_JUNCTION_ANCHOR` (8): an `N` is trusted as a real junction only with >= 8 matched bases each side. Higher trusts
  fewer junctions; lower lets coincidental short anchors count as splices.
- `SPLICE_FLANKING_DELETION_MAX_BP` (5): a `D <= 5` straddling an exon boundary is absorbed into the `N`. Higher folds
  bigger boundary deletions into the intron; lower keeps them as an explicit `D`.
- `DEFAULT_MAX_SHIFT` (5): max bp an intron slides to reach a higher splice-motif tier (canonicalize). Higher repositions
  junctions further onto a canonical motif (risk: a coincidental motif); lower leaves more junctions off-motif.

Confidence and noise floors (constants):

- `RESCUE_MAPQ` (60): MAPQ given on a swap or a rescued unique placement. The BAM unique-mapper convention; past 60 is
  out of spec, lower makes downstream MAPQ filters discount correctly-rescued reads.
- `SUPP_RESCUE_MAPQ_CAP` (55): ceiling on a primary+supp merge, flagging a constructed alignment. Higher toward 60 makes
  merges look bwa-native; lower discounts real rescued splices more.
- `PRIMARY_AS_UNMAP_THRESHOLD` (30): a primary still below this AS after liftback is unmapped. Higher unmaps more
  borderline reads; lower (toward 19) keeps the short-anchor noise bwa's `-T 19` surfaced.
- `SUPP_AS_DROP_THRESHOLD` (30): a supplementary still below this AS is dropped. Same trade as the primary floor.

Upstream bwa-mem2 flags (not tars config, but liftback depends on them):

- `-T 19` (score floor): set below bwa's default 30 so bwa emits the short-anchor supps rescue merges. Raising it to 30
  starves rescue; tars re-imposes the 30 floor after rescue via the AS floors above.
- `-h 75` (XA cap): a read over the cap is emitted MAPQ 0 with no XA; tars unmaps such a genomic primary
  (`exceedsMappingCap`). Higher keeps more high-multiplicity reads as listed multimappers; lower unmaps more.

Scoring: bwa-mem affine scores used to re-score lifted boundaries are `MATCH / MISMATCH = +1 / -4` and
`GAP_OPEN / GAP_EXTEND = -6 / -1` (soft-clips not penalised). `NH = max(numLoci, 1)` counts distinct genomic loci, not
emitted records.

</details>

<details>
<summary><strong>Edge cases</strong></summary>

These are the rare exits and refinements. On the whole-sample run the unmap/lift-fail exits below total 22,899
records, a tiny fraction, but each prevents a specific class of wrong placement.

**Over-XA-cap unmap (too many genomic loci).** bwa is run with `-h 75` (XA cap). A read mapping to more loci than the cap
is emitted at MAPQ 0 with **no** XA (the alt list is suppressed, not truncated). The over-cap rule unmaps such a primary:
too many genomic places to trust.

```
condition = (inputMapq == 0 AND numXaAlts == 0 AND comp == REF_ONLY)
```

The `REF_ONLY` gate is load-bearing. A tx-contig primary can hit 75+ transcript contigs that all lift to **one** genomic
locus, so its suppressed XA must not be read as too many genomic places. Only a genomic (REF_ONLY) primary is unmapped;
TX_ONLY / REF_AND_TX lift normally. Keyed on the input state, not the post-lift MAPQ: with no XA the resolver would
otherwise see a single locus and rescue MAPQ to 60.

**Excluded-region contamination (-rna_unmap_regions).** Curated rRNA / 7SL / acrocentric / multi-map contamination
zones. A read lifting into one is contamination: the primary is **unmapped** (kept, not dropped), a supplementary is
**dropped**. Checked **post-lift** on genomic coordinates, not pre-lift. A tx-contig read's input coords are `chrN_tx`
and cannot be tested against the genomic region list, and some excluded zones (acrocentric p-arms) are themselves in the
transcriptome, so contamination reaches them via tx contigs and is only visible once lifted. The primary is unmapped by
flipping its result to UNMAPPED (so the mate is coordinated via the cache); dropped supps have their SA entry removed
from the primary's SA tag.

**AS floor unmap / drop (residual short-anchor noise).** bwa runs at `-T 19` to surface short-anchor supps for rescue.
After liftback, records still below the default `-T 30` AS floor that were not rescued/extended/collapsed are residual
noise.

* `PRIMARY_AS_UNMAP_THRESHOLD = 30`: a primary still below 30 is **unmapped** (not dropped, to keep the pair + SA
intact). Gated on rescue running; AS is never recomputed, so post-processed primaries keep a stale-low AS and are exempt.
* `SUPP_AS_DROP_THRESHOLD = 30`: a supplementary in the [19, 30) AS band that survived rescue + tail-extend is **dropped**.

**Flanking-deletion absorption in translation.** A small `D` straddling an exon boundary is usually a tx-FASTA off-by-N
artefact. `ContigTranslator` folds a `D <= SPLICE_FLANKING_DELETION_MAX_BP (5)` into the intron `N`; a larger `D` (6+ bp)
is preserved after the `N`. Bound is inclusive.

```
90M 3D 61M  with the 3D straddling a boundary  ->  90M 5000N 61M    (3D absorbed)
90M 8D 61M                                     ->  90M 5000N 8D 61M (8D preserved)
```

**LIFT_FAILED and anchor floors.**

* A position outside any transcript span (inter-transcript spacer, or past the contig end) -> LIFT_FAILED: read marked
unmapped, ref/start/CIGAR/MAPQ cleared, XA and SA stripped.
* A supplementary whose own lift failed is instead mirrored onto its primary's lifted coords, keeping the 0x800 flag.
* Terminal micro-anchors below the floors are folded into a softclip during translation only when adjacent to an
existing softclip: bare-anchor floor 1 bp, softclip-adjacent floor 3 bp. A flush terminal anchor (read ends on the
anchor, no softclip) is left for the Step 1 collapse pass, which scores it against the genome. (Both are distinct
from `MIN_JUNCTION_ANCHOR` = 8, the discriminator/collapse trust floor.)
* Read overhang past the first/last exon span is clamped to a leading/trailing softclip.

**Hidden-tie MAPQ block.** The single-locus MAPQ-0-to-60 rescue does **not** fire when there is an unresolved hidden
tie: `XS == AS` on a ref-only primary outside any annotated exon. That signals a genuine ambiguous placement, so the
MAPQ stays at 0.

</details>
