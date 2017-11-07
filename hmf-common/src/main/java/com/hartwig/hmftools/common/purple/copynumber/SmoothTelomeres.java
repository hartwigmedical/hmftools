package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;

import org.jetbrains.annotations.NotNull;

class SmoothTelomeres {

    @NotNull
    static List<FittedRegion> smooth(@NotNull final List<FittedRegion> regions) {
        if (regions.size() > 1) {
            final FittedRegion telomere = regions.get(0);
            final FittedRegion neighbour = regions.get(1);
            if (Doubles.isZero(telomere.tumorCopyNumber())) {
                final List<FittedRegion> result = Lists.newArrayList();
                final FittedRegion fixedTelomere =
                        ImmutableFittedRegion.builder().from(telomere).tumorCopyNumber(neighbour.tumorCopyNumber()).build();
                result.add(fixedTelomere);
                result.addAll(regions.subList(1, regions.size()));
                return result;
            }
        }

        return regions;
    }
}
