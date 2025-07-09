package com.hartwig.hmftools.geneutils.paneldesign;

public class PanelBuilderConstants
{
    // Common constants.
    public static final int PROBE_LENGTH = 120;
    // Above this quality score we may consider a candidate probe to be acceptable without checking other candidates.
    public static final double PROBE_QUALITY_ACCEPT = 0.5;
    // All probes with quality score below this value are unconditionally rejected.
    public static final double PROBE_QUALITY_REJECT = 0.08;
    // When covering a region with probes, how many bases in the region are allowed to not be covered?
    public static final int UNCOVERED_BASES_MAX = 10;

    // Target genes constants.
    // Space between the gene and upstream/downstream region probes.
    public static final int GENE_UPDOWNSTREAM_GAP = 1000;
    // Region within which the upstream/downstream probes are selected.
    public static final int GENE_UPDOWNSTREAM_REGION = 2000;
    // Only probe intronic regions when the gene has <= this number of exons.
    public static final int GENE_MAX_EXONS_TO_ADD_INTRON = 19;
    // Space between exon and intron probe.
    public static final int GENE_EXON_FLANK_GAP = 1000;
    // Intron probes are selected within a region of this many bases.
    public static final int GENE_EXON_FLANK_REGION = 1000;
    // Introns shorter than this many bases are ignored.
    public static final int GENE_MIN_INTRON_LENGTH = GENE_EXON_FLANK_GAP * 2 + GENE_EXON_FLANK_REGION;
    // Introns shorter than this many bases get 1 probe, longer introns get 2 probes.
    public static final int GENE_LONG_INTRON_LENGTH = GENE_EXON_FLANK_GAP * 2 + GENE_EXON_FLANK_REGION * 3;

    // Copy number backbone constants.
    public static final int CN_BACKBONE_PARTITION_SIZE = 1_000_000;
    public static final int CN_BACKBONE_CENTROMERE_MARGIN = 5_000_000;
    // Aiming to pick heterozygous sites which are common in the population.
    public static final double CN_BACKBONE_GNOMAD_FREQ_MIN = 0.3;
    public static final double CN_BACKBONE_GNOMAD_FREQ_MAX = 0.7;
    // Very high quality threshold because there are many sites to pick from to we may as well get the best.
    public static final double CN_BACKBONE_QUALITY_MIN = 0.8;

    // GC content bounds for probes used for determining copy number.
    // These are tight bounds because different GC can affect the probe amplification process which will affect the calculated copy number.
    public static final double CN_GC_TARGET = 0.45;
    public static final double CN_GC_TOLERANCE = 0.025;
    public static final double CN_GC_MIN = CN_GC_TARGET - CN_GC_TOLERANCE;
    public static final double CN_GC_MAX = CN_GC_TARGET + CN_GC_TARGET;

    // Output naming.
    public static final String PANEL_PROBES_FILE = "panel_probes.tsv";
    public static final String REJECTED_REGIONS_FILE = "rejected_regions.tsv";
}
