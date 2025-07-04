package com.hartwig.hmftools.common.mappability;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.BaseRegion;

public class ProbeQualityWindow extends BaseRegion
{
    private final float mQualityScore;

    public ProbeQualityWindow(int posStart, int posEnd, float qualityScore)
    {
        super(posStart, posEnd);
        this.mQualityScore = qualityScore;
    }

    public float getQualityScore()
    {
        return mQualityScore;
    }

    @Override
    public String toString()
    {
        return format("%d-%d quality(%f)", start(), end(), mQualityScore);
    }

    @Override
    public boolean equals(final Object obj)
    {
        return super.equals(obj) && ((ProbeQualityWindow) obj).mQualityScore == mQualityScore;
    }
}
