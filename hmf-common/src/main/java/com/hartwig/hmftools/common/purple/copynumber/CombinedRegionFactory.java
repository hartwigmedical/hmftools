package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ModifiableFittedRegion;

import org.jetbrains.annotations.NotNull;

class CombinedRegionFactory {

    @NotNull
    static CombinedRegion reduce(@NotNull final CombinedRegion parent, long start, long end) {
        assert (start >= parent.start());
        assert (end <= parent.end());

        int bafCount = 0;
        int observedTumorRatioCount = 0;
        for (FittedRegion fittedRegion : parent.regions()) {
            if (fittedRegion.start() >= start && fittedRegion.end() <= end) {
                bafCount += fittedRegion.bafCount();
                observedTumorRatioCount += fittedRegion.observedTumorRatioCount();
            }
        }

        final ModifiableFittedRegion smallerRegion = ModifiableFittedRegion.create()
                .from(parent.region())
                .setStart(start)
                .setEnd(end)
                .setBafCount(bafCount)
                .setObservedTumorRatioCount(observedTumorRatioCount);

        CombinedRegion result = new CombinedRegion(parent.isBafWeighted(), smallerRegion, false);
        result.setCopyNumberMethod(parent.copyNumberMethod());

        for (FittedRegion fittedRegion : parent.regions()) {
            if (fittedRegion.start() >= start && fittedRegion.end() <= end) {
                result.extend(fittedRegion);
            }
        }

        return result;
    }

}
