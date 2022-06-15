package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;

final class BreakendSelector {

    private BreakendSelector() {
    }

    @NotNull
    public static List<LinxBreakend> selectInterestingUnreportedBreakends(@NotNull List<LinxBreakend> allBreakends,
            @NotNull List<LinxFusion> reportableFusions, @NotNull KnownFusionCache knownFusionCache) {
        List<LinxBreakend> interestingUnreportedBreakends = Lists.newArrayList();
        for (LinxBreakend breakend : allBreakends) {
            if (!breakend.reportedDisruption() && breakend.disruptive()) {
                boolean isBreakInPromiscuousExonRange = isBreakInPromiscuousExonRange(knownFusionCache, breakend);

                if (isBreakInPromiscuousExonRange) {
                    interestingUnreportedBreakends.add(breakend);
                }
            }
        }
        return interestingUnreportedBreakends;
    }

    private static boolean isBreakInPromiscuousExonRange(@NotNull KnownFusionCache knownFusionCache, @NotNull LinxBreakend breakend) {
//        knownFusionCache.getDataByType(KnownFusionType.PROMISCUOUS_3)
        return false;
    }
}
