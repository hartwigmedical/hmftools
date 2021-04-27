package com.hartwig.hmftools.purple.region;

import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;

public interface FittedRegionFactory
{
    @NotNull
    FittedRegion fitRegion(final double purity, final double normFactor, final @NotNull ObservedRegion observedRegion);

    @NotNull
    List<FittedRegion> fitRegion(final double purity, final double normFactor, @NotNull final Collection<ObservedRegion> observedRegions);
}
