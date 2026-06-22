# TARS

**TARS** (Transcript Alignment for RNA Splicing) produces splice-aware RNA alignments using `bwa-mem2`. RNA reads
are aligned against the genome FASTA with annotated multi-exon transcript contigs (`*_tx`) appended. Liftback then
rewrites those transcript-contig alignments back to genomic coordinates, so the output is an ordinary genomic RNA BAM:
no `*_tx` in the header or on records, XA/SA/mate fields genomic, and spliced reads carried as `N` CIGAR ops.

The lifted BAM feeds into REDUX for normal (splice-agnostic) dedup, then ISOFOX.

## Contents

* [The RNA flow](#the-rna-flow)
* [Commands](#commands)
  + [Building transcript contigs](#building-transcript-contigs-splicefastabuilder)
  + [Running liftback](#running-liftback-spliceliftback)
  + [Output](#output)
* [Flags](#flags)
* [How a read gets lifted](#how-a-read-gets-lifted)
  + [Step 0: translate the transcript CIGAR to genome](#step-0-translate-the-transcript-cigar-to-genome)
  + [Step 1: ref vs tx discriminator](#step-1-ref-vs-tx-discriminator)
  + [Step 2: rescue via supplementary](#step-2-rescue-via-supplementary)
  + [Step 3: terminal micro-junction collapse](#step-3-terminal-micro-junction-collapse)
  + [Step 4: tail extension](#step-4-tail-extension)
  + [Step 5: canonicalize](#step-5-canonicalize)
* [Modifications done to each record](#modifications-done-to-each-record)
* [Edge cases](#edge-cases)
* [Constants reference](#constants-reference)
* [Tests](#tests)

## The RNA flow

```
  ensembl  ->  SpliceFastaBuilder  ->  genome.fasta + *_tx contigs  ->  bwa-mem2 index
  (offline, once per release)          + sidecar TSV

  RNA fastq  ->  bwa-mem2 align (name-grouped, no sort)  ->  SpliceLiftBack  ->  genomic RNA BAM (sorted) ->  REDUX  ->  ISOFOX
```

1. **`SpliceFastaBuilder`** (offline, once per ensembl release) builds the transcript-contig FASTA and a sidecar TSV
   mapping packed transcript intervals to gene/transcript/exon spans.
2. Concatenate the transcript contigs onto the genome FASTA and re-index with `bwa-mem2`.
3. Align RNA reads with `bwa-mem2` against the combined FASTA. bwa's default output is **name-grouped** (each
   fragment's mates and supplementaries are contiguous, in FASTQ order). Feed it to liftback directly, no sort needed.
4. **`SpliceLiftBack`** consumes that name-grouped BAM in a single pass (no input sort, no index, no fragment cache),
   lifts each fragment to genomic coordinates, and writes a coord-sorted and indexed genomic BAM.
5. Feed the lifted BAM into REDUX for normal dedup.

## Commands

### Building transcript contigs (SpliceFastaBuilder)

```
java -cp tars.jar com.hartwig.hmftools.tars.fasta.SpliceFastaBuilder
    -ensembl_data_dir /ref_data/ensembl_data_cache/38/
    -ref_genome /path_to_fasta/genome.fasta
    -ref_genome_version V38
    -output_dir /path_to_output/
```

Outputs `ref_genome_v38_rna_contigs.fasta` and `ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv`. Concatenate the
FASTA onto the genome FASTA and `bwa-mem2 index` the result before aligning.

### Running liftback (SpliceLiftBack)

```
java -cp tars.jar com.hartwig.hmftools.tars.liftback.SpliceLiftBack
    -input_bam SAMPLE_ID.bwa_tx.namegrouped.bam
    -ref_genome /path_to_fasta/genome_plus_tx.fasta
    -ref_genome_version V38
    -contig_sidecar /path_to/ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv
    -rna_unmap_regions /ref_data/rna/38/rna_excluded_regions.38.tsv
    -bamtool /path_to_samtools/
    -output_dir /path_to_output/
    -output_id SAMPLE_ID
    -threads 24
```

### Output

`<output_id>.bam` (+ `.bai`), coord-sorted and ready for REDUX. The summary always writes to `*.liftback.summary.tsv`;
the per-record `*.liftback.records.tsv` / `*.liftback.alignments.tsv` are written only with `-write_liftback_tsv`.

## Flags

**Required**

| Flag               | Description                                                                  |
|--------------------|------------------------------------------------------------------------------|
| input_bam          | bwa-mem2 output against the combined FASTA, **name-grouped** (not coord-sorted)|
| ref_genome         | The same combined genome + transcript FASTA used at alignment                 |
| ref_genome_version | `V37` or `V38`                                                                |
| contig_sidecar     | Contig sidecar TSV from `SpliceFastaBuilder` (`*.rna_contigs_mappings.tsv`)    |
| bamtool            | samtools / sambamba path, used to decompress the input and sort + index output|
| output_dir         | Directory for the lifted BAM and summary                                      |

**Optional**

| Flag               | Default | Description                                                              |
|--------------------|---------|--------------------------------------------------------------------------|
| output_id          | (none)  | Output filename prefix; also namespaces per-worker shard intermediates   |
| rna_unmap_regions  | (none)  | Curated regions (rRNA / 7SL / multi-map); a read lifting into one is unmapped|
| write_liftback_tsv | off     | Emit per-record debug TSVs; off by default (whole-sample TSVs run to 100s of GB)|
| threads            | 1       | Lift-worker count; the lift phase scales with this (the final sort tail does not). Run many (prod uses 24)|
| log_level          | INFO    | `log_debug` / `log_level` for more detail                                |

**Tuning thresholds**

Junction rescue and softclip tail-extension both **run by default** (no on/off flag); these expose their thresholds for
tuning the false-positive / yield tradeoff. Defaults are also in the [constants reference](#constants-reference).

| Flag                        | Default   | Description                                                       |
|-----------------------------|-----------|------------------------------------------------------------------|
| rescue_min_intron           | 21        | Min intron length for a primary+supp merge                       |
| rescue_max_intron           | 1000000   | Max intron length for a merge                                    |
| rescue_min_anchor_overhang  | 3         | Min matched bases each side of a merged junction                 |
| rescue_max_chain_depth      | 4         | Max supplementary pieces merged into one primary                 |
| rescue_softclip_tolerance   | 5         | Primary/supp overlap tolerance when snapping to a junction       |
| rescue_max_boundary_shift   | 8         | Max over-extended bases trimmed when probing a boundary          |
| rescue_min_partial_match_run| 11        | Min exon-proximal matched run for a partial ref-verify rescue    |
| tail_min_softclip           | 3         | Min terminal softclip length to consider                         |
| tail_min_extension          | 3         | Min ref-matching bases converted to M                            |
| tail_max_extension          | 30        | Max bases walked into a softclip                                 |

Note: no `ensembl_data_dir` - liftback reads exon/junction annotation from the sidecar (only `SpliceFastaBuilder` needs ensembl).

## How a read gets lifted

Liftback translates each tx-contig alignment to genomic coordinates, chooses the primary between the ref and tx
placements, then refines the spliced CIGAR through four passes (rescue, collapse, tail-extend, canonicalize), each
feeding the next. The CIGAR examples below use a 5000 bp intron for illustration; real intron lengths vary.

### Step 0: translate the transcript CIGAR to genome

`ContigTranslator` walks the transcript-contig CIGAR through the transcript's exons and inserts an `N` at each exon
boundary the read crosses. This is the core splice transform: a read aligned flat against a `*_tx` contig becomes a
spliced genomic read.

```
  tx contig:   chr7_ENST..._tx   151M        (flat against the spliced transcript, no gaps)
  exon boundary falls at read base 90, exons are 5000 bp apart on the genome
  ->  genome:  chr7  90M 5000N 61M
```

A read that sits entirely inside one exon crosses no boundary, so it lifts with its CIGAR unchanged (plain `M`, no `N`).

### Step 1: ref vs tx discriminator

bwa-mem2 can align a read to a ref contig, a transcript (`tx`) contig, or both, at one or several genomic loci.
`LiftBackDiscriminator` categorizes that alignment set and picks the primary alignment (and which XA alts are worth
keeping). Two definitions drive every rule:

- **tx has a real junction** = at least one `N` whose both flanking `M` runs are >= 8 bp (`MIN_JUNCTION_ANCHOR`); one
  strong junction trusts the whole read, so a short anchor elsewhere does not disqualify a multi-exon read.
- **ref full match** = the ref (genomic) alt has no softclip, i.e. bwa aligns the whole read contiguously at that locus
  (mismatches/SNVs allowed - NM is not checked). The genome explaining the read with no junction points to it being
  genomic there (retained intron, intronic, or a single-exon transcript) rather than spliced.

#### No-swap cases (93% of primaries)

`REF_SINGLE` + `REF_MULTI` (71.5%) - ref only, no tx match (one locus or multimapped): no change. Ref alts are already
genomic, so they pass through with original coords and MAPQ; the discriminator does nothing.
```
  151M @ chr1   (single)   or   151M @ chr1 + 151M @ chr9   (multi)   ->  unchanged
```
`TX_SINGLE` (12.3%) - tx only, one locus (the core splice case):
```
  90M 5000N 61M   ->  lifts as the spliced read
```
`BOTH_AGREE` (9.3%) - ref and tx, same CIGAR, no N:
```
  ref 151M  +  tx 151M   ->  agree; keep, drop the duplicate tx alt
```

#### Swaps: correcting bwa's primary pick (2.31% of primaries)

In the four `BOTH_*` contest categories the read has a ref placement and a tx placement at one locus (or paralogous
loci) and bwa chose the wrong one as primary. The discriminator swaps it: the winner's coordinates move into the
record's main fields, and the displaced placement is either kept as an alternate in the read's `XA` tag (SAM's
secondary-alignment list, `XA:Z:chrom,pos,CIGAR,NM;...`) or dropped. bwa picked wrong for a different reason each time.

**tx wins** - the read is genuinely spliced, so the tx placement becomes primary and the displaced ref is kept in `XA`.

`BOTH_TX_JUNCTION_REF_SOFTCLIP` (0.58%) - tx spans a real junction; ref softclipped at the boundary, never crossed it.
bwa made the softclipped ref primary because it favours a placement near the read's mate (a proper pair) over a
higher-scoring tx hit on a different contig.
```
  ref alt:   45S 106M          (clipped at the boundary)
  tx alt:    45M 5000N 106M    (real junction, both anchors >= 8 bp)
  ->  primary becomes the tx alt, MAPQ 60; old ref demoted to XA
```
`BOTH_MULTI_TX_JUNCTION` (1.69%) - tx has a real junction; the ref alts are intronless paralogs at other loci. bwa
picked a paralog because the intronless copy aligns contiguously and scores as high as the spliced parent.
```
  tx alt:    90M 5000N 61M  @ chr1    (spliced parent gene)
  ref alt:   151M           @ chr15   (intronless paralog)
  ->  primary becomes the tx alt; the paralog ref demoted to XA
```

**ref wins** - the read is unspliced, so the tx N-CIGAR is the artifact (retained intron / pre-mRNA / DNA); ref becomes
primary and the spurious tx alt is dropped.

`BOTH_TX_JUNCTION_REF_MATCH` (<0.01%) - ref reads cleanly straight through the supposed intron. bwa made the tx alt
primary (a clean match against the spliced contig), but the genome explains the read with no gap, so ref is right.
```
  ref alt:   151M              (no softclip)
  tx alt:    90M 5000N 61M
  ->  primary becomes the ref alt; tx alt dropped
```
`BOTH_TX_SOFTCLIP_REF_MATCH` (0.04%) - tx only softclips at the boundary (never spanned a junction); ref reads through.
The read sits inside a retained intron, which only the genome explains.
```
  tx alt:    90M 61S           (softclips where the exon ends, no N)
  ref alt:   151M              (no softclip, contiguous)
  ->  primary becomes the ref alt; tx alt dropped
```

#### Reference

Full category distribution (exp9 whole sample, 442.9M primaries):

| Category                       |     Records | % primaries | Outcome                                       |
|--------------------------------|------------:|------------:|-----------------------------------------------|
| REF_SINGLE                     | 251,677,114 |   56.82%    | genomic placement, lifts as-is                |
| REF_MULTI                      |  64,827,831 |   14.64%    | multimapped ref, no swap                      |
| TX_SINGLE                      |  54,662,124 |   12.34%    | core splice case, tx-only -> `M N M`          |
| BOTH_AGREE                     |  41,366,113 |    9.34%    | ref and tx agree, no contest                  |
| BOTH_MULTI                     |  12,132,293 |    2.74%    | genuine multimapper, no swap                  |
| BOTH_MULTI_TX_JUNCTION         |   7,504,615 |    1.69%    | multi-locus, **tx wins** (ref alts intronless)|
| TX_MULTI                       |   5,308,785 |    1.20%    | multimapped tx, no swap                       |
| BOTH_AMBIGUOUS                 |   2,698,912 |    0.61%    | no rule fires, no swap                        |
| BOTH_TX_JUNCTION_REF_SOFTCLIP  |   2,559,367 |    0.58%    | single-locus contest, **tx wins**             |
| BOTH_TX_SOFTCLIP_REF_MATCH     |     175,313 |    0.04%    | single-locus contest, **ref wins**            |
| BOTH_TX_JUNCTION_REF_MATCH     |       5,616 |   <0.01%    | single-locus contest, **ref wins**            |

The no-swap categories: `REF_SINGLE` / `REF_MULTI` (ref only), `TX_SINGLE` (tx only), `BOTH_AGREE` (identical), and the
three problem cases that keep bwa's pick - `BOTH_MULTI` (genuine multimapper), `TX_MULTI` (multimapped tx), and
`BOTH_AMBIGUOUS` (no rule fires).

Full decision tree:
```
loci >= 2 ?
|
+-- YES (multi-locus)
|     hasRef AND hasTx ?
|     |  YES:  TxHasNCigar AND NOT RefHasNCigar -> BOTH_MULTI_TX_JUNCTION   [tx faithful; ref alts intronless]
|     |        else                             -> BOTH_MULTI               [genuine multimapper, no swap]
|     |  tx only                                -> TX_MULTI
|     |  ref only                               -> REF_MULTI
|
+-- NO (single locus)
      only ref                                  -> REF_SINGLE
      only tx                                   -> TX_SINGLE
      both ref and tx:
         same CIGAR and no N                    -> BOTH_AGREE                    [no contest]
         TxHasNCigar AND RefSoftClipped         -> BOTH_TX_JUNCTION_REF_SOFTCLIP *** TX WINS ***
         TxHasNCigar AND RefFullMatch           -> BOTH_TX_JUNCTION_REF_MATCH    *** REF WINS ***
         TxSoftClipAtBoundary AND RefFullMatch  -> BOTH_TX_SOFTCLIP_REF_MATCH    *** REF WINS ***
         otherwise                              -> BOTH_AMBIGUOUS                [no rule fires, no swap]

  TxHasNCigar   = >=1 N with both flanks >= 8 bp        RefFullMatch = ref alt with no softclip
```

### Step 2: rescue via supplementary

Runs by default. bwa runs at a low score floor (`-T 19`) on purpose, which surfaces short-anchor supplementary
alignments. Rescue merges a primary's terminal softclip with such a supplementary across an annotated (or
motif-supported) junction into one `M N M` primary, and drops the merged supplementary.

```
  primary:   106M 45S         softclip and...
  supp:      45M               ...a short anchor across the junction
  ->  merged primary:  106M 5000N 45M    (supplementary dropped, MAPQ capped at 55)
```

Gated on softclip complementarity, intron length [21 .. 1,000,000], anchor overhang >= 3 bp, coverage overlap, and
chain depth <= 4. The snap point is chosen by splice-motif tier.

### Step 3: terminal micro-junction collapse

`TerminalMicroJunctionCollapser` removes a fabricated tiny terminal anchor (< 8 bp) sitting across an intron when the
contiguous genome explains the tail at least as well as bwa's split.

```
  in:   146M 2000N 5M         (5 bp terminal anchor across a "junction")
  contiguous genome matches straight through:
  ->    151M                  (the anchor was not a real junction)
```

The decision is a head-to-head bwa-mem score: anchor-on-the-far-exon vs the whole terminal window extended
contiguously. The score-maximising prefix is reclaimed into the near exon as `M`, the rest stays softclipped. A real
short junction (anchor scores strictly higher on the far exon than the contiguous walk) is kept.

### Step 4: tail extension

Runs by default. Walk a terminal softclip into contiguous genome with no `N` (intron retention). Softclip must be
>= 3 bp, extension lands in [3 .. 30] bp.

```
  in:   131M 20S              (terminal clip that actually matches the genome)
  ->    151M                  (fully reclaimed)   or   139M 12S  (partial, 12 bp residual stays clipped)
```

Skipped when an annotated intron starts inside the window and the walk stalls short: that range belongs to rescue.

### Step 5: canonicalize

`JunctionCanonicalizer` slides an intron up to 5 bp onto a higher splice-motif tier
(NONE < SEMI_CANONICAL < CANONICAL < ANNOTATED), choosing the smallest shift that strictly improves the tier with every
moved base still matching. CIGAR only: the intron position moves but the alignment start never does.

```
  in:   90M 5000N 61M         (junction 2 bp off a GT-AG motif)
  ->    92M 4998N 59M         (slid 2 bp onto the canonical motif)
```

## Modifications done to each record

On a successful lift (`LiftBackRecordOps` / `MateFieldPatcher`):

1. Reference name, alignment start, CIGAR, strand and MAPQ are set to the resolved genomic placement.
2. XA rebuilt to genomic alt coordinates; SA rewritten to genomic; MC set to the mate's lifted CIGAR.
3. `XS:A:+/-` written only when the lifted CIGAR has an `N` **and** the transcript strand is known.
4. MD always dropped (never rebuilt); NM recomputed against the genomic reference.
5. `NH = max(numLoci, 1)`, counting distinct genomic loci among best-scoring alignments only, so one junction across
   many transcript contigs does not inflate it.
6. Mate fields patched: signed 5'-to-5' TLEN, mate coords, unmapped-mate parking.

**MAPQ** is carried from bwa and never raised by collapse / tail-extend / canonicalize. It is set to 60 when the
discriminator swaps to tx, or when a single-locus read had MAPQ 0 with no unresolved hidden tie. A primary+supp rescue
merge caps it at 55.

## Edge cases

These are the rare exits and refinements. On the exp9 whole-sample run the unmap/lift-fail exits below total 20,009
records, **0.004%** of all, but each prevents a specific class of wrong placement.

<details>
<summary>Over-XA-cap unmap (too many genomic loci)</summary>

bwa is run with `-h 75` (XA cap). A read mapping to more loci than the cap is emitted at MAPQ 0 with **no** XA (the alt
list is suppressed, not truncated). The over-cap rule unmaps such a primary: too many genomic places to trust.

```
condition = (inputMapq == 0 AND numXaAlts == 0 AND comp == REF_ONLY)
```

The `REF_ONLY` gate is load-bearing. A tx-contig primary hits 75+ transcript contigs of one gene (a shared exon) that
all lift to **one** genomic locus, so its suppressed XA must not be read as too many genomic places. Only a genomic
(REF_ONLY) primary is unmapped; TX_ONLY / REF_AND_TX lift normally. Keyed on the input state, not the post-lift MAPQ:
with no XA the resolver would otherwise see a single locus and rescue MAPQ to 60.
</details>

<details>
<summary>Excluded-region contamination (-rna_unmap_regions)</summary>

Curated rRNA / 7SL / acrocentric / multi-map contamination zones. A read lifting into one is contamination: the primary
is **unmapped** (REDUX-style, kept not dropped), a supplementary is **dropped**.

Checked **post-lift** on genomic coordinates, not pre-lift. A tx-contig read's input coords are `chrN_tx` and cannot be
tested against the genomic region list, and some excluded zones (acrocentric p-arms) are themselves in the
transcriptome, so contamination reaches them via tx contigs and is only visible once lifted. The primary is unmapped by
flipping its result to UNMAPPED (so the mate is coordinated via the cache); dropped supps have their SA entry removed
from the primary's SA tag.
</details>

<details>
<summary>AS floor unmap / drop (residual short-anchor noise)</summary>

bwa runs at `-T 19` to surface short-anchor supps for rescue. After liftback, records still below the default `-T 30`
AS floor that were not rescued/extended/collapsed are residual noise.

* `PRIMARY_AS_UNMAP_THRESHOLD = 30`: a primary still below 30 is **unmapped** (not dropped, to keep the pair + SA
  intact). Gated on rescue running; AS is never recomputed, so post-processed primaries keep a stale-low AS and are exempt.
* `SUPP_AS_DROP_THRESHOLD = 30`: a supplementary in the [19, 30) AS band that survived rescue + tail-extend is **dropped**.
</details>

<details>
<summary>Flanking-deletion absorption in translation</summary>

A small `D` straddling an exon boundary is usually a tx-FASTA off-by-N artefact. `ContigTranslator` folds a
`D <= SPLICE_FLANKING_DELETION_MAX_BP (5)` into the intron `N`; a larger `D` (6+ bp) is preserved after the `N`. Bound
is inclusive.

```
  90M 3D 61M  with the 3D straddling a boundary  ->  90M 5000N 61M    (3D absorbed)
  90M 8D 61M                                     ->  90M 5000N 8D 61M (8D preserved)
```
</details>

<details>
<summary>LIFT_FAILED and anchor floors</summary>

* A position outside any transcript span (inter-transcript spacer, or past the contig end) -> LIFT_FAILED: read marked
  unmapped, ref/start/CIGAR/MAPQ cleared, XA and SA stripped.
* A supplementary whose own lift failed is instead mirrored onto its primary's lifted coords, keeping the 0x800 flag.
* Terminal micro-anchors below the floors are folded into a softclip during translation: bare-anchor floor 1 bp,
  softclip-adjacent floor 3 bp. (Distinct from `MIN_JUNCTION_ANCHOR = 8`, the discriminator/collapse trust floor.)
* Read overhang past the first/last exon span is clamped to a leading/trailing softclip.
</details>

<details>
<summary>Hidden-tie MAPQ block</summary>

The single-locus MAPQ-0-to-60 rescue does **not** fire when there is an unresolved hidden tie: `XS == AS` on a ref-only
primary outside any annotated exon. That signals a genuine ambiguous placement, so the MAPQ stays at 0.
</details>

## Constants reference

The important ones (all in `TarsConstants`; the rescue/tail values are also tunable via the [flags](#flags) above):

| Constant                        | Value     | What it gates                                                      |
|---------------------------------|-----------|-------------------------------------------------------------------|
| `MIN_JUNCTION_ANCHOR`           | 8         | an `N` is a trusted junction only if at least one has both flanks >= 8 bp |
| `SPLICE_FLANKING_DELETION_MAX_BP`| 5        | a `D <= 5` across an exon boundary is absorbed into the `N`        |
| `RESCUE_MAPQ`                   | 60        | MAPQ assigned on a discriminator swap or a rescued unique placement |
| `SUPP_RESCUE_MAPQ_CAP`          | 55        | ceiling on a primary+supp rescue merge (flags a constructed alignment) |
| `PRIMARY_AS_UNMAP_THRESHOLD`    | 30        | primary still below this AS after liftback is unmapped             |
| `SUPP_AS_DROP_THRESHOLD`        | 30        | supplementary still below this AS is dropped                       |
| `DEFAULT_MIN_INTRON_LENGTH`     | 21        | shortest intron a primary+supp rescue may bridge                  |
| `DEFAULT_MAX_INTRON_LENGTH`     | 1,000,000 | longest intron a rescue may bridge                                |
| `DEFAULT_MIN_ANCHOR_OVERHANG`   | 3         | matched bases needed each side of a merged junction               |
| `DEFAULT_MAX_CHAIN_DEPTH`       | 4         | supplementary pieces merged into one primary                      |
| `DEFAULT_MAX_EXTENSION`         | 30        | max bases a terminal softclip is walked into contiguous genome     |
| `DEFAULT_MAX_SHIFT`             | 5         | max bp an intron slides to reach a higher splice-motif tier        |

bwa-mem affine scoring used to re-score lifted boundaries: `MATCH / MISMATCH = +1 / -4`, `GAP_OPEN / GAP_EXTEND = -6 / -1`
(soft-clips not penalised). `NH = max(numLoci, 1)` (counts distinct genomic loci, not emitted records).

## Tests

End-to-end behaviour is pinned by `LiftBackEndToEndTest` (full per-group engine):

| Scenario                                       | What it pins                                                  |
|------------------------------------------------|--------------------------------------------------------------|
| exonSpanningReadLiftsToJunctionCigar           | tx read spanning two exons -> `M N M`, mate fields patched   |
| unliftableReadIsMarkedUnmapped                 | position past the contig end -> unmapped, mate flagged       |
| supplementarySaTagRewrittenToGenomicCoords     | split read; supp + SA lifted to genomic coords               |
| terminalSoftclipExtendedIntoContiguousGenome   | intron-retention tail walked into the genome (tail-extend)   |
| splitReadRescuedAcrossAnnotatedJunction        | primary softclip + short-anchor supp merged into one `M N M` |
| offMotifJunctionSlidToCanonical                | junction slid 2 bp onto a GT-AG motif (canonicalize)         |
| nativeGenomicReadPassesThroughUnchanged        | a read already on the genome is left untouched               |

Component layers below the engine: `LiftBackResolverTest` (categorize + resolved result), `SpliceLiftBackApplyTest`
(record mutation: CIGAR / XA / XS / NH / unmap), `LiftBackDiscriminatorTest` (full category matrix + swap/drop), and
the rescue / tail-extend suites (individual refinement passes).
</content>
</invoke>
