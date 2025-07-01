package com.hartwig.hmftools.common.mappability;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.BaseRegion;

public class ProbeQualityWindow extends BaseRegion
{
    private final float mQualityScore;

    public ProbeQualityWindow(final int posStart, final int posEnd, final float qualityScore)
    {
        super(posStart, posEnd);
        this.mQualityScore = qualityScore;
    }

    public float getQualityScore()
    {
        return mQualityScore;
    }

    @Override
    public String toString() { return format("%d-%d quality(%f)", start(), end(), mQualityScore); }
}
