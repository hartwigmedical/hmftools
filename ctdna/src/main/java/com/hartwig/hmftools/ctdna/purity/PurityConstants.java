package com.hartwig.hmftools.ctdna.purity;

public class PurityConstants
{
    public static final double MAX_SUBCLONAL_LIKELIHOOD = 0.1;
    public static final double MAX_REPEAT_COUNT = 3;

    public static final int MIN_QUAL_PER_AD = 18;
    public static final int LOW_QUAL_NOISE_CUTOFF = 25;
    public static final double SAMPLE_ALLELE_OUTLIER_PROBABIILTY = 0.00001;

    public static final double MAX_COPY_NUMBER = 6;
    public static final double CLONAL_COPY_NUMBER_MARGIN = 0.2;

    public static final double DEFAULT_NOISE_READS_PER_MILLION = 30;
    public static final double DEFAULT_NOISE_READS_PER_MILLION_DUAL_STRAND = 1;

    public static final String PURPLE_CTDNA_SOMATIC_VCF_ID = ".purple.somatic.ctdna.";
}
