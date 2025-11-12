package com.hartwig.hmftools.common.utils.pcf;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class CobaltSegment extends ChrBaseRegion
{
    public final double MeanRatio;

    public CobaltSegment(final String chromosome, final int posStart, final int posEnd, final double meanRatio)
    {
        super(chromosome, posStart, posEnd);
        MeanRatio = meanRatio;
    }
}
