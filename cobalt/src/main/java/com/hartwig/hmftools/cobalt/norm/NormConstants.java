package com.hartwig.hmftools.cobalt.norm;

import static com.hartwig.hmftools.common.genome.gc.GCProfile.MIN_MAPPABLE_PERCENTAGE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.WINDOW_SIZE;

public class NormConstants
{
    public static final int REGION_SIZE = WINDOW_SIZE;

    public static final double MAPPABILITY_THRESHOLD = MIN_MAPPABLE_PERCENTAGE;

    public static final double MIN_ENRICHMENT_RATIO = 0.1;
}
