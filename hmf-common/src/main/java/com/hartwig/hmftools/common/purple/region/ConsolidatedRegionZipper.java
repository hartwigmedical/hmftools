package com.hartwig.hmftools.common.purple.region;

import java.util.List;
import java.util.function.BiFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedRegion;
import com.hartwig.hmftools.common.purple.ImmutableFittedRegion;
import com.hartwig.hmftools.common.zipper.RegionZipper;
import com.hartwig.hmftools.common.zipper.RegionZipperHandler;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ConsolidatedRegionZipper implements RegionZipperHandler<ConsolidatedRegion, FittedRegion> {

    @NotNull
    private final BiFunction<ConsolidatedRegion, FittedRegion, FittedRegion> transform;
    private final List<FittedRegion> result = Lists.newArrayList();

    @Nullable
    private ConsolidatedRegion consolidatedRegion;

    private ConsolidatedRegionZipper(
            @NotNull final BiFunction<ConsolidatedRegion, FittedRegion, FittedRegion> transform) {
        this.transform = transform;
    }

    @NotNull
    private List<FittedRegion> result() {
        return result;
    }

    @Override
    public void left(final ConsolidatedRegion region) {
        consolidatedRegion = region;
    }

    @Override
    public void right(final FittedRegion copyNumber) {
        if (consolidatedRegion != null && consolidatedRegion.overlaps(copyNumber)) {
            result.add(transform.apply(consolidatedRegion, copyNumber));
        } else {
            result.add(copyNumber);
        }
    }

    @NotNull
    public static List<FittedRegion> insertSmoothRegions(@NotNull final List<ConsolidatedRegion> smoothRegions,
            @NotNull final List<FittedRegion> fittedRegions) {
        BiFunction<ConsolidatedRegion, FittedRegion, FittedRegion> transform = (consolidatedRegion, copyNumber) -> ImmutableFittedRegion
                .builder().from(
                copyNumber).segmentBAF(consolidatedRegion.averageObservedBAF()).segmentTumorCopyNumber(
                consolidatedRegion.averageTumorCopyNumber()).build();

        ConsolidatedRegionZipper myThing = new ConsolidatedRegionZipper(transform);
        RegionZipper.zip(smoothRegions, fittedRegions, myThing);
        return myThing.result();
    }

    @NotNull
    public static List<FittedRegion> insertHighConfidenceRegions(
            @NotNull final List<ConsolidatedRegion> highConfidenceRegions, @NotNull final List<FittedRegion> fittedRegions) {
        BiFunction<ConsolidatedRegion, FittedRegion, FittedRegion> transform = (consolidatedRegion, copyNumber) -> ImmutableFittedRegion
                .builder().from(
                copyNumber).broadBAF(consolidatedRegion.averageObservedBAF()).broadTumorCopyNumber(
                consolidatedRegion.averageTumorCopyNumber()).build();

        ConsolidatedRegionZipper myThing = new ConsolidatedRegionZipper(transform);
        RegionZipper.zip(highConfidenceRegions, fittedRegions, myThing);
        return myThing.result();
    }
}
