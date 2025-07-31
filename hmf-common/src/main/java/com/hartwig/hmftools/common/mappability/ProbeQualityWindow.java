package com.hartwig.hmftools.common.mappability;

import com.hartwig.hmftools.common.region.BaseRegion;

public record ProbeQualityWindow(
        BaseRegion region,
        float qualityScore
) implements Comparable<ProbeQualityWindow>
{
    public ProbeQualityWindow
    {
        if(!(qualityScore >= 0 && qualityScore <= 1))
        {
            throw new IllegalArgumentException("qualityScore must be in the range [0, 1]");
        }
    }

    public ProbeQualityWindow(int posStart, int posEnd, float qualityScore)
    {
        this(new BaseRegion(posStart, posEnd), qualityScore);
    }

    @Override
    public int compareTo(final ProbeQualityWindow other)
    {
        return region.compareTo(other.region);
    }
}
