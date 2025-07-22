package com.hartwig.hmftools.geneutils.paneldesign;

public class PanelBuilderConstants
{
    // Common constants.
    public static final int PROBE_LENGTH = 120;
    // By default, have wide GC content tolerance since it matters less for general probes.
    public static final double GENERAL_GC_TARGET = 0.45;
    public static final double GENERAL_GC_TOLERANCE = 1;    // I.e. any GC is ok
    // When covering a large region with probes, how many bases a probe may be shifted left or right from its "ideal" tiled position.
    public static final int PROBE_SHIFT_MAX = 5;

    // Target genes constants.
    // All the region and gap sizes should be much larger than the probe size to avoid probe overlap.
    // Space between the gene and upstream/downstream region probes.
    public static final int GENE_UPDOWNSTREAM_GAP = 1000;
    // Region within which the upstream/downstream probes are selected.
    public static final int GENE_UPDOWNSTREAM_REGION = 2000;
    // Coding exons are expanded by this many bases each side to include splice points in probe coverage.
    public static final int GENE_CODING_REGION_EXPAND = 10;
    public static final double GENE_EXON_QUALITY_MIN = 0.05;
    // Minimum space between exon and intron probe.
    public static final int GENE_EXON_FLANK_GAP = 1000;
    // Exon flank probes are selected within a region of this many bases.
    public static final int GENE_EXON_FLANK_REGION_MIN = 1000;
    public static final int GENE_EXON_FLANK_REGION_MAX = 5000;
    public static final double GENE_CN_QUALITY_MIN = 0.5;

    // Copy number backbone constants.
    public static final int CN_BACKBONE_PARTITION_SIZE = 1_000_000;
    public static final int CN_BACKBONE_CENTROMERE_MARGIN = 5_000_000;
    // Aiming to pick heterozygous sites which are common in the population.
    public static final double CN_BACKBONE_GNOMAD_FREQ_MIN = 0.3;
    public static final double CN_BACKBONE_GNOMAD_FREQ_MAX = 0.7;
    // Very high probe quality threshold because there are many sites to pick from to we may as well get the best.
    public static final double CN_BACKBONE_QUALITY_MIN = 0.8;

    // GC content bounds for probes used for determining copy number.
    // These are tight bounds because different GC can affect the probe amplification process which will affect the calculated copy number.
    public static final double CN_GC_TARGET = 0.45;
    public static final double CN_GC_TOLERANCE = 0.05;
    // Early stopping on the probe search if the GC is within this tolerance. (Performance tuning parameter.)
    public static final double CN_GC_OPTIMAL_TOLERANCE = 0.001;

    // Custom regions parameters.
    public static final double CUSTOM_REGION_QUALITY_MIN = 0.1;

    // Output naming.
    public static final String PANEL_PROBES_FILE_STEM = "panel_probes";
    public static final String TARGET_REGIONS_FILE_NAME = "target_regions.bed";
    public static final String REJECTED_REGIONS_FILE_STEM = "rejected_regions";
    // This output can get very large (multiple GB) so write it in compressed format.
    public static final String CANDIDATE_PROBES_FILE_NAME = "candidate_probes.tsv.gz";

    static
    {
        if(!(PROBE_LENGTH >= 50))
        {
            // Wouldn't recommend going much smaller than this because it was designed to work with 120b probe
            throw new IllegalArgumentException();
        }
        if(!(PROBE_SHIFT_MAX >= 0 && PROBE_SHIFT_MAX < PROBE_LENGTH))
        {
            throw new IllegalArgumentException();
        }
        if(!(GENE_UPDOWNSTREAM_REGION > PROBE_LENGTH))
        {
            throw new IllegalArgumentException();
        }
        if(!(GENE_EXON_FLANK_GAP > GENE_CODING_REGION_EXPAND))
        {
            throw new IllegalArgumentException();
        }
        if(!(GENE_EXON_FLANK_REGION_MIN > PROBE_LENGTH))
        {
            throw new IllegalArgumentException();
        }
        if(!(GENE_EXON_FLANK_REGION_MAX > GENE_EXON_FLANK_REGION_MIN))
        {
            throw new IllegalArgumentException();
        }
        if(!(CN_BACKBONE_PARTITION_SIZE >= PROBE_LENGTH * 2))
        {
            throw new IllegalArgumentException();
        }
    }
}
