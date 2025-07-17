package com.hartwig.hmftools.common.mappability;

import com.hartwig.hmftools.common.region.BaseRegion;

public record ProbeQualityWindow(
        BaseRegion region,
        float qualityScore
) implements Comparable<ProbeQualityWindow>
{
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
