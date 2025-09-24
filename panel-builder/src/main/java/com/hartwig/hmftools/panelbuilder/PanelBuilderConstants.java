package com.hartwig.hmftools.panelbuilder;

public class PanelBuilderConstants
{
    public static final String APP_NAME = "PanelBuilder";

    // Common constants.
    public static final int PROBE_LENGTH = 120;
    // If we can't calculate the quality score of a probe, what quality score do we assign?
    public static final double DEFAULT_PROBE_QUALITY = 0;

    // Parameters when covering an arbitrary region with probes.
    // How many bases a probe may be shifted left or right from its "ideal" tiled position.
    public static final int PROBE_SHIFT_MAX = 5;
    // How many bases may be uncovered before we add another probe?
    // Uncovered bases will always be on the edges of the regions.
    // This is desirable because reads from sequencing will include an area slightly larger than the probe.
    public static final int REGION_UNCOVERED_MAX = 20;
    // When probes cover more than the target region, how much to allocate to probe overlap vs. extension outside the target region.
    // 0 = maximise overlap. 1 = maximise extension.
    public static final double PROBE_OVERLAP_EXTENSION_BALANCE = 0.5;

    // Gene probes constants.
    // All the region and gap sizes should be much larger than the probe size to avoid probe overlap.
    // Space between the gene and upstream/downstream region probes.
    public static final int GENE_UPDOWNSTREAM_GAP = 1000;
    // Region within which the upstream/downstream probes are selected.
    public static final int GENE_UPDOWNSTREAM_REGION = 2000;
    // Coding exons are expanded by this many bases each side to include splice points in probe coverage.
    public static final int GENE_CODING_REGION_EXPAND = 10;
    // Minimum space between exon and intron probe.
    public static final int GENE_EXON_FLANK_GAP = 1000;
    // Exon flank probes are selected within a region of this many bases.
    public static final int GENE_EXON_FLANK_REGION_MIN = 1000;
    public static final int GENE_EXON_FLANK_REGION_MAX = 5000;
    // How many bases upstream to cover as part of the promoter region.
    public static final int GENE_PROMOTER_REGION = 500;
    public static final double GENE_GENERAL_QUALITY_MIN = 0.05;
    public static final double GENE_GENERAL_GC_TARGET = 0.45;
    public static final double GENE_GENERAL_GC_TOLERANCE = 1;
    public static final double GENE_CN_QUALITY_MIN = 0.5;

    // Copy number backbone constants.
    public static final int CN_BACKBONE_RESOLUTION_KB_DEFAULT = 1_000;
    // Region excluded, per side of the centromere.
    public static final int CN_BACKBONE_CENTROMERE_MARGIN = 3_000_000;
    // Aiming to pick heterozygous sites which are common in the population.
    public static final double CN_BACKBONE_GNOMAD_FREQ_MIN = 0.3;
    public static final double CN_BACKBONE_GNOMAD_FREQ_MAX = 0.7;
    // Very high probe quality threshold because there are many sites to pick from to we may as well get the best.
    public static final double CN_BACKBONE_QUALITY_MIN = 0.8;
    public static final double CN_BACKBONE_ALTERNATE_QUALITY_MIN = 1.0;

    // CDR3 probe constants.
    // Very low quality score threshold because many of these CDR3 regions are similar to another, which is ok.
    // We manually vetted the region list to ensure the "off-target" is other CDR3 regions.
    public static final double CDR3_QUALITY_MIN = 0.01;
    public static final double CDR3_GC_TARGET = 0.45;
    public static final double CDR3_GC_TOLERANCE = 1;

    // GC content bounds for probes used for determining copy number.
    // These are tight bounds because different GC can affect the probe amplification process which will affect the calculated copy number.
    public static final double CN_GC_TARGET = 0.45;
    public static final double CN_GC_TOLERANCE = 0.05;
    // Early stopping on the probe search if the GC is within this tolerance. (Performance tuning parameter.)
    public static final double CN_GC_OPTIMAL_TOLERANCE = 0.005;

    // Custom regions parameters.
    public static final double CUSTOM_REGION_QUALITY_MIN = 0.1;
    public static final double CUSTOM_REGION_GC_TARGET = 0.45;
    public static final double CUSTOM_REGION_GC_TOLERANCE = 1;

    // Sample variants parameters.
    public static final int SAMPLE_PROBES_MAX_DEFAULT = 500;
    public static final double SAMPLE_NONDRIVER_QUALITY_MIN = 0.1;
    public static final double SAMPLE_NONDRIVER_GC_TARGET = 0.45;
    public static final double SAMPLE_NONDRIVER_GC_TOLERANCE = 0.15;
    // Driver criteria is less strict because it's more important to include them compared to nondrivers.
    public static final double SAMPLE_DRIVER_QUALITY_MIN = 0.05;
    public static final double SAMPLE_DRIVER_GC_TARGET = 0.45;
    public static final double SAMPLE_DRIVER_GC_TOLERANCE = 1;
    public static final double SAMPLE_VAF_MIN = 0.05;
    public static final int SAMPLE_FRAGMENT_COUNT_MIN = 11;
    public static final int SAMPLE_REPEAT_COUNT_MAX = 3;
    public static final int SAMPLE_SV_BREAKENDS_PER_GENE_MAX = 5;
    public static final int SAMPLE_INDEL_LENGTH_MAX = 31;
    public static final double SAMPLE_SUBCLONAL_LIKELIHOOD_MIN = 0.95;
    // If a variant probe differs from the ref genome by more than this many bases, consider it a novel sequence, and keep the probe even if
    // it overlaps with existing probes.
    public static final int SAMPLE_NOVEL_SEQUENCE_BASES_MIN = 5;

    // Output naming.
    public static final String PANEL_PROBES_FILE_STEM = "panel_probes";
    public static final String TARGET_REGIONS_FILE_NAME = "target_regions.bed";
    public static final String REJECTED_REGIONS_FILE_STEM = "rejected_regions";
    public static final String CANDIDATE_REGIONS_FILE_NAME = "candidate_regions.bed";
    // This output can get very large (multiple GB) so write it in compressed format.
    public static final String CANDIDATE_PROBES_FILE_NAME = "candidate_probes.tsv.gz";
    public static final String GENE_STATS_FILE_NAME = "gene_stats.tsv";

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
        if(!(0 <= PROBE_OVERLAP_EXTENSION_BALANCE && PROBE_OVERLAP_EXTENSION_BALANCE <= 1))
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
        if(!(GENE_UPDOWNSTREAM_GAP > GENE_PROMOTER_REGION))
        {
            throw new IllegalArgumentException();
        }
    }
}
