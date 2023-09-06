package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class BreakendUtil
{
    @NotNull
    public static List<Pair<LinxBreakend, LinxBreakend>> createPairsPerSvId(@NotNull List<LinxBreakend> breakends)
    {
        Map<Integer, Pair<LinxBreakend, LinxBreakend>> mapPerSvId = Maps.newHashMap();
        for(LinxBreakend breakend : breakends)
        {
            if(!mapPerSvId.containsKey(breakend.svId()))
            {
                LinxBreakend paired = findPaired(breakends, breakend);
                if(paired != null)
                {
                    mapPerSvId.put(breakend.svId(), Pair.of(breakend, paired));
                }
            }
        }
        return Lists.newArrayList(mapPerSvId.values());
    }

    @Nullable
    private static LinxBreakend findPaired(@NotNull List<LinxBreakend> breakends, @NotNull LinxBreakend breakendToFindPairFor)
    {
        for(LinxBreakend breakend : breakends)
        {
            if(breakend != breakendToFindPairFor && breakend.svId() == breakendToFindPairFor.svId())
            {
                return breakend;
            }
        }
        return null;
    }
}
