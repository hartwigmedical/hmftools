package com.hartwig.hmftools.common.purple.copynumber;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ModifiableFittedRegion;

import org.jetbrains.annotations.NotNull;

class ExtractChildren {

    @NotNull
    static List<CombinedRegion> fillGaps(@NotNull final CombinedRegion parent, final List<CombinedRegion> children) {
        if (children.isEmpty()) {
            return Collections.singletonList(parent);
        }

        final List<CombinedRegion> result = Lists.newArrayList();
        long nextStart = parent.start();
        for (CombinedRegion child : children) {
            if (child.start() > nextStart) {
                result.add(reduce(parent, nextStart, child.start() - 1));
            }

            result.add(child);
            nextStart = child.end() + 1;
        }

        if (nextStart <= parent.end()) {
            result.add(reduce(parent, nextStart, parent.end()));
        }

        return result;
    }

    @NotNull
    private static CombinedRegion reduce(@NotNull final CombinedRegion parent, long start, long end) {
        assert (start >= parent.start());
        assert (end <= parent.end());

        int bafCount = 0;
        int depthWindowCount = 0;
        for (FittedRegion fittedRegion : parent.regions()) {
            if (fittedRegion.start() >= start && fittedRegion.end() <= end) {
                bafCount += fittedRegion.bafCount();
                depthWindowCount += fittedRegion.depthWindowCount();
            }
        }

        final ModifiableFittedRegion smallerRegion = ModifiableFittedRegion.create()
                .from(parent.region())
                .setStart(start)
                .setEnd(end)
                .setBafCount(bafCount)
                .setDepthWindowCount(depthWindowCount);

        CombinedRegion result = new CombinedRegionImpl(smallerRegion);
        result.setCopyNumberMethod(parent.copyNumberMethod());

        for (FittedRegion fittedRegion : parent.regions()) {
            if (fittedRegion.start() >= start && fittedRegion.end() <= end) {
                result.extend(fittedRegion);
            }
        }

        return result;
    }

}
