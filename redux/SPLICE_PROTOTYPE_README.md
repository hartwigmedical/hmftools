# Splice-aware alignment prototype

Experimental helpers under `com.hartwig.hmftools.redux.splice` that make `bwa-mem2` splice-aware in an effort to replace STAR.

The flow is:

1. **SpliceFastaBuilder**: offline, run once per ensembl release. Emits a FASTA of synthetic per-transcript contigs (concatenated exon sequences) plus a sidecar TSV mapping each contig back to gene / transcript / chromosome / exon spans.
2. Append the synthetic contigs to the genome FASTA and re-index with `bwa-mem2`.
3. Align reads with `bwa-mem2` against the combined FASTA.
4. **SpliceLiftBack**: per sample. Reads the sidecar, translates transcript-coordinate alignments back to genome coordinates (CIGAR with `N` gaps across introns).

## Tool 1: SpliceFastaBuilder (offline, one-time per ensembl)

```
java -cp redux/target/redux-2.0-jar-with-dependencies.jar \
    com.hartwig.hmftools.redux.splice.SpliceFastaBuilder \
    -ensembl_data_dir ~/Documents/hartwig/common-resources-public/ensembl_data_cache/38 \
    -ref_genome ~/Documents/hartwig/data/rna/ref_genome/38/GRCh38_masked_exclusions_alts_hlas.fasta \
    -ref_genome_version V38 \
    -output_dir ~/Documents/hartwig/temp-output
```

Outputs `splice.transcript_contigs.fasta` + `splice.transcript_contigs.sidecar.tsv` in the output dir.
