package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.abs;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class MatchedSineElement
{
    public final String RmId;
    public final BaseRegion StartRegion;
    public final BaseRegion EndRegion;

    public final List<SvBreakendData> Breakends;

    public MatchedSineElement(final String rmId, final BaseRegion startRegion, final BaseRegion endRegion)
    {
        RmId = rmId;
        StartRegion = startRegion;
        EndRegion = endRegion;
        Breakends = Lists.newArrayList();
    }

    public boolean isProximatePosition(final int position, final int maxDistance)
    {
        if(abs(StartRegion.start() - position) <= maxDistance)
            return true;

        if(abs(StartRegion.end() - position) <= maxDistance)
            return true;

        if(abs(EndRegion.start() - position) <= maxDistance)
            return true;

        if(abs(EndRegion.end() - position) <= maxDistance)
            return true;

        return false;
    }
}
