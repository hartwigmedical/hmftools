package com.hartwig.hmftools.geneutils.paneldesign;

public class PanelConstants
{
    // copy number backbone
    public static final int CN_BACKBONE_PARTITION_SIZE = 1_000_000;
    public static final int CN_BACKBONE_CENTROMERE_MARGIN = 5_000_000;
    public static final double CN_BACKBONE_MAPPABILITY = 1;
    public static final double CN_BACKBONE_GNMOD_FREQ_MIN = 0.3;
    public static final double CN_BACKBONE_GNMOD_FREQ_MAX = 0.7;
    public static final double CN_BACKBONE_GC_RATIO_MIN = 0.45;

    // common
    public static final int PROBE_LENGTH = 120;

    // target genes
    public static final double PROBE_GC_MIN = 0.35;
    public static final double PROBE_GC_MAX = 0.55;
    public static final int GENE_MAX_CANDIDATE_PROBES = 8;
    public static final int GENE_MIN_INTRON_LENGTH = 3000;
    public static final int GENE_LONG_INTRON_LENGTH = 5000;
    public static final int GENE_MAX_EXONS_TO_ADD_INTRON = 19;
    public static final int GENE_FLANKING_DISTANCE = 1000;
    public static final int GENE_CANDIDATE_REGION_SIZE = PROBE_LENGTH * GENE_MAX_CANDIDATE_PROBES;

    // BlastN settings
    public static final int BLASTN_WORD_SIZE = 15;
    public static final int MIN_BLAST_ALIGNMENT_LENGTH = 30;

    // BlastN results
    public static final double MAX_PROBE_SUM_BLASTN_BITSCORE = 2500;

}
