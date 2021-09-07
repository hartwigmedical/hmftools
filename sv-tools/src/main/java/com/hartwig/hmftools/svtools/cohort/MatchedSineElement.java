package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;

import java.util.List;

import com.google.common.collect.Lists;

public class MatchedSineElement
{
    public final RepeatMaskerData RmStart;
    public final RepeatMaskerData RmEnd;

    public final List<SvBreakendData> Breakends;

    public MatchedSineElement(final RepeatMaskerData rmStart, final RepeatMaskerData rmEnd)
    {
        RmStart = rmStart;
        RmEnd = rmEnd;
        Breakends = Lists.newArrayList();
    }

    public String combinedRmId() { return String.format("%s_%s", RmStart.RmId, RmEnd.RmId); }

    public boolean isPositionWithin(final int position)
    {
        return positionWithin(position, RmStart.Region.start(), RmEnd.Region.end());
    }
}
