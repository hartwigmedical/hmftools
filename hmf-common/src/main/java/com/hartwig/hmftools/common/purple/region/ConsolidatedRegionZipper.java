package com.hartwig.hmftools.common.purple.region;

import java.util.List;
import java.util.function.BiFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedCopyNumber;
import com.hartwig.hmftools.common.purple.ImmutableFittedCopyNumber;
import com.hartwig.hmftools.common.zipper.RegionZipper;
import com.hartwig.hmftools.common.zipper.RegionZipperHandler;

import org.jetbrains.annotations.Nullable;

public class ConsolidatedRegionZipper implements RegionZipperHandler<ConsolidatedRegion, FittedCopyNumber> {

    private final BiFunction<ConsolidatedRegion, FittedCopyNumber, FittedCopyNumber> transform;
    private final List<FittedCopyNumber> result = Lists.newArrayList();

    @Nullable
    private ConsolidatedRegion consolidatedRegion;

    private ConsolidatedRegionZipper(
            final BiFunction<ConsolidatedRegion, FittedCopyNumber, FittedCopyNumber> transform) {
        this.transform = transform;
    }

    private List<FittedCopyNumber> getResult() {
        return result;
    }

    @Override
    public void left(final ConsolidatedRegion region) {
        consolidatedRegion = region;
    }

    @Override
    public void right(final FittedCopyNumber copyNumber) {
        if (consolidatedRegion != null && consolidatedRegion.overlaps(copyNumber)) {
            result.add(transform.apply(consolidatedRegion, copyNumber));
        } else {
            result.add(copyNumber);
        }

    }

    public static List<FittedCopyNumber> insertSmoothRegions(List<ConsolidatedRegion> megaRegions,
            List<FittedCopyNumber> copyNumbers) {

        BiFunction<ConsolidatedRegion, FittedCopyNumber, FittedCopyNumber> transform = (consolidatedRegion, copyNumber) -> ImmutableFittedCopyNumber
                .builder()
                .from(copyNumber)
                .segmentBAF(consolidatedRegion.averageObservedBAF())
                .segmentTumorCopyNumber(consolidatedRegion.averageTumorCopyNumber())
                .build();

        ConsolidatedRegionZipper myThing = new ConsolidatedRegionZipper(transform);
        RegionZipper.zip(megaRegions, copyNumbers, myThing);
        return myThing.getResult();
    }

    public static List<FittedCopyNumber> insertHighConfidenceRegions(List<ConsolidatedRegion> megaRegions,
            List<FittedCopyNumber> copyNumbers) {

        BiFunction<ConsolidatedRegion, FittedCopyNumber, FittedCopyNumber> transform = (consolidatedRegion, copyNumber) -> ImmutableFittedCopyNumber
                .builder()
                .from(copyNumber)
                .broadBAF(consolidatedRegion.averageObservedBAF())
                .broadTumorCopyNumber(consolidatedRegion.averageTumorCopyNumber())
                .build();

        ConsolidatedRegionZipper myThing = new ConsolidatedRegionZipper(transform);
        RegionZipper.zip(megaRegions, copyNumbers, myThing);
        return myThing.getResult();
    }

}
