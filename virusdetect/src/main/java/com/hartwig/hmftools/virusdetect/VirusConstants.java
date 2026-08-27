package com.hartwig.hmftools.virusdetect;

public class VirusConstants
{
    public static final String APP_NAME = "VirusDetect";

    // Provisional defaults, to be tuned on real data.

    // Minimum soft-clip length for a read to count as a candidate.
    public static final int MIN_SOFT_CLIP_BASES_DEFAULT = 20;

    // Minimum BWA alignment score, as a fraction of read length, to keep an alignment.
    public static final double MIN_ALIGNMENT_SCORE_FRACTION_DEFAULT = 0.5;

    // Candidate reads submitted to BWA per alignment call (bounds heap).
    public static final int ALIGNMENT_BATCH_SIZE_DEFAULT = 100_000;

    // Minimum reads for an oncology group to be considered present. Not user-exposed yet; tuned internally first.
    public static final int MIN_READS_DEFAULT = 2;
}
