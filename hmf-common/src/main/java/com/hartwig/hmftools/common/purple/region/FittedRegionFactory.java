package com.hartwig.hmftools.common.purple.region;

import java.util.Collection;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public interface FittedRegionFactory {

    @NotNull
    List<FittedRegion> fitRegion(final double purity, final double normFactor,
            @NotNull final Collection<ObservedRegion> observedRegions);

    @NotNull
    FittedRegion fitRegion(final double purity, final double normFactor, final @NotNull ObservedRegion observedRegion);

}
