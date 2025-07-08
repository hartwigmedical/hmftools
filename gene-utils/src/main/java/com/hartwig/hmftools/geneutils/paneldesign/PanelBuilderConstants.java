package com.hartwig.hmftools.geneutils.paneldesign;

public class PanelBuilderConstants
{
    // Common
    public static final int PROBE_LENGTH = 120;
    public static final double PROBE_GC_MIN = 0.35;
    public static final double PROBE_GC_MAX = 0.55;
    public static final double PROBE_QUALITY_SCORE_MIN = 0.3;

    // Target genes
    public static final int GENE_MAX_CANDIDATE_PROBES = 8;
    public static final int GENE_MIN_INTRON_LENGTH = 3000;
    public static final int GENE_LONG_INTRON_LENGTH = 5000;
    public static final int GENE_MAX_EXONS_TO_ADD_INTRON = 19;
    public static final int GENE_FLANKING_DISTANCE = 1000;
    public static final int GENE_CANDIDATE_REGION_SIZE = PROBE_LENGTH * GENE_MAX_CANDIDATE_PROBES;

    // Copy number backbone
    public static final int CN_BACKBONE_PARTITION_SIZE = 1_000_000;
    public static final int CN_BACKBONE_CENTROMERE_MARGIN = 5_000_000;
    public static final double CN_BACKBONE_MAPPABILITY_MIN = 1;
    public static final double CN_BACKBONE_GNOMAD_FREQ_MIN = 0.3;
    public static final double CN_BACKBONE_GNOMAD_FREQ_MAX = 0.7;
    public static final double CN_BACKBONE_GC_RATIO_MIN = 0.45;

    // Output naming
    public static final String PANEL_PROBES_FILE = "panel_probes.tsv";
    public static final String REJECTED_REGIONS_FILE = "rejected_regions.tsv";
}
