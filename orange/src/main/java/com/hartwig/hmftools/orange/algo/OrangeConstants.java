package com.hartwig.hmftools.orange.algo;

import static com.hartwig.hmftools.datamodel.driver.DriverInterpretation.DRIVER_LIKELIHOOD_LOW_THRESHOLD;

import com.hartwig.hmftools.common.purple.GermlineStatus;

public final class OrangeConstants
{
    // Purple
    public static final double PURPLE_ARM_CN_GAIN_THRESHOLD = 1.25;
    public static final double PURPLE_ARM_CN_LOSS_THRESHOLD = 0.75;

    public static final String PURPLE_ARM_CN_LOSS = "LOSS";
    public static final String PURPLE_ARM_CN_GAIN = "GAIN";
    public static final String PURPLE_ARM_CN_DIPLOID = GermlineStatus.DIPLOID.toString();

    public static final String PURPLE_AMP_DEL_PARTIAL = "PARTIAL";
    public static final String PURPLE_AMP_DEL_FULL = "FULL";

    // RNA
    public static final double HIGH_EXPRESSION_PERCENTILE_CUTOFF = 0.99;
    public static final double LOW_EXPRESSION_PERCENTILE_CUTOFF = 0.01;
    public static final int MAX_EXPRESSION_GENE_COUNT = 100; // TEMP

    public static boolean isCandidateLikelihood(final double likelihood) { return likelihood <= DRIVER_LIKELIHOOD_LOW_THRESHOLD; }
}
