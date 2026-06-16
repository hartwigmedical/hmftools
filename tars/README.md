# TARS

**TARS** (Transcript Alignment for RNA Splicing) produces splice-aware RNA alignments using bwa-mem2. RNA reads
are aligned against the genome FASTA with annotated multi-exon transcript contigs (`*_tx`) appended; liftback rewrites
those tx-contig alignments back to genomic coordinates so the output is an ordinary genomic RNA BAM: no `*_tx` in the
header or on records, XA/SA/mate fields genomic, spliced reads carried as `N` CIGAR ops.

TARS is a standalone tool (not a REDUX dedup stage). The end-to-end RNA flow is:

1. **`SpliceFastaBuilder`** (offline, once per ensembl release) builds the transcript-contig FASTA and a sidecar TSV
   mapping packed transcript intervals to gene/transcript/exon spans.
2. Append the transcript contigs to the genome FASTA and re-index with `bwa-mem2`.
3. Align RNA reads with `bwa-mem2` against the combined FASTA. bwa's default output is **name-grouped** (each
   fragment's mates + supplementaries contiguous, FASTQ order) - feed it to liftback directly, no sort needed.
4. **`SpliceLiftBack`** consumes that name-grouped BAM in a single pass (no input sort, no index, no fragment cache),
   lifts each fragment to genomic coordinates, and writes a coord-sorted + indexed genomic BAM.
5. Feed the lifted BAM into REDUX for normal dedup. REDUX itself is splice-agnostic.

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

### Command

```
java -cp tars.jar com.hartwig.hmftools.tars.liftback.SpliceLiftBack 
    -input_bam SAMPLE_ID.bwa_tx.namegrouped.bam 
    -ref_genome /path_to_fasta/genome_plus_tx.fasta 
    -ref_genome_version V38 
    -contig_sidecar /path_to/ref_genome_v38_rna_contigs.rna_contigs_mappings.tsv 
    -ensembl_data_dir /ref_data/ensembl_data_cache/38/ 
    -rna_unmap_regions /ref_data/rna/38/rna_excluded_regions.38.tsv 
    -rescue_via_supp 
    -extend_softclip_tails 
    -bamtool /path_to_samtools/ 
    -output_dir /path_to_output/ 
    -output_id SAMPLE_ID 
    -threads 24
```

| Argument              | Required | Description                                                                      |
|-----------------------|----------|----------------------------------------------------------------------------------|
| input_bam             | Required | bwa-mem2 output against the combined FASTA, **name-grouped** (not coord-sorted)  |
| ref_genome            | Required | The same combined genome + transcript FASTA used at alignment                    |
| contig_sidecar        | Required | Contig sidecar TSV from `SpliceFastaBuilder` (`*.rna_contigs_mappings.tsv`)       |
| ensembl_data_dir      | Required | Ensembl cache: annotated exons + junctions, used by the discriminator and rescue |
| bamtool               | Required | samtools / sambamba path, used to decompress the input and sort + index output   |
| rna_unmap_regions     | Optional | Curated regions (rRNA / 7SL / multi-map); a fragment is dropped before lifting   |
| rescue_via_supp       | Optional | Merge primary + supplementary pieces into one spliced primary                    |
| extend_softclip_tails | Optional | Recover terminal softclipped bases that match the reference (intron retention)   |
| write_liftback_tsv    | Optional | Emit per-record debug TSVs; off by default (whole-sample TSVs run to 100s of GB) |

Output is `<output_id>.bam` (+ `.bai`), coord-sorted and ready for REDUX. The summary always writes to
`*.liftback.summary.tsv`; the per-record `*.liftback.records.tsv` / `*.liftback.alignments.tsv` are written only
with `-write_liftback_tsv`.

### Translation

`ContigTranslator` walks a tx-contig CIGAR through the transcript's exons, inserting an `N` at each exon boundary
crossed. Small `D` ops next to a new `N` are absorbed into it; terminal anchors under 3 bp are folded into softclip.

### Ref-vs-tx discriminator

When a read has both a ref alignment and a tx (N-carrying) alignment at one locus, `LiftBackDiscriminator` picks the
primary. An `N` is a real junction only if both flanking M runs are >= 8 bp. Ref wins when it matches cleanly across
the supposed intron (retained intron / pre-mRNA / DNA contamination), or tx only softclips at the boundary. Tx wins
when it has a real junction and ref is softclipped (ref never spanned the junction).

### Post-process passes

Run per mate-group in this load-bearing order, each feeding the next:

1. **Rescue** (`JunctionRescueResolver`) — merge a primary's terminal softclip with a short-anchor supplementary across
   a junction into one `M N M` primary. Gated on softclip complementarity, intron length, anchor overhangs and coverage
   overlap; snap point chosen by motif tier. Merged supps are dropped; the primary's MAPQ is capped at 55.
2. **Collapse** (`TerminalMicroJunctionCollapser`) — drop a fabricated tiny terminal anchor (<= 3 bp) across an intron.
   The anchor bases must match the contiguous genome exactly; reclaimed softclip bases past the anchor tolerate
   max(1, length/10) mismatches. A real short anchor has an intron between, so it doesn't match contiguously and is kept.
3. **Extend** (`SoftclipTailExtender`) — walk a terminal softclip into contiguous genome (intron retention, no `N`),
   min 3 bp, capped at 30 bp.
4. **Canonicalize** (`JunctionCanonicalizer`) — slide an intron up to 5 bp onto a higher splice motif (GT-AG / CT-AC),
   smallest shift that strictly improves the tier with every moved base still matching. CIGAR only; start never moves.

### MAPQ, NH and mutation scope

MAPQ is carried from bwa and never raised on collapse/extend/canonicalize. It is set to 60 when the discriminator swaps
to tx, or a single-locus read had MAPQ 0 with no unresolved hidden tie; rescue caps it at 55. NH = `max(numLoci, 1)`,
counting distinct genomic loci among best-scoring alignments only, so one junction across many transcript contigs does
not inflate it. All record mutation lives in `LiftBackRecordOps`, `MateFieldPatcher` and `LiftBackGroupProcessor`;
everything else is pure. MD is always dropped (not rebuilt), NM is recomputed against the genomic reference, and read
bases are never reverse-complemented when the strand flag flips (bwa already oriented the read).
