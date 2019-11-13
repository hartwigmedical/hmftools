package com.hartwig.hmftools.purple;

import java.util.List;
import java.util.function.BiFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.utils.zipper.RegionZipper;
import com.hartwig.hmftools.common.utils.zipper.RegionZipperHandler;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class PurpleRegionZipper implements RegionZipperHandler<PurpleCopyNumber, FittedRegion> {

    @NotNull
    private final BiFunction<PurpleCopyNumber, FittedRegion, FittedRegion> transform;
    private final List<FittedRegion> result = Lists.newArrayList();

    @Nullable
    private PurpleCopyNumber consolidatedRegion;

    private PurpleRegionZipper(@NotNull final BiFunction<PurpleCopyNumber, FittedRegion, FittedRegion> transform) {
        this.transform = transform;
    }

    @NotNull
    private List<FittedRegion> result() {
        return result;
    }

    @Override
    public void enterChromosome(@NotNull final String chromosome) {
        // Empty
    }

    @Override
    public void primary(@NotNull final PurpleCopyNumber region) {
        consolidatedRegion = region;
    }

    @Override
    public void secondary(@NotNull final FittedRegion copyNumber) {
        if (consolidatedRegion != null && consolidatedRegion.overlaps(copyNumber)) {
            result.add(transform.apply(consolidatedRegion, copyNumber));
        } else {
            result.add(copyNumber);
        }
    }

    static List<FittedRegion> updateRegionsWithCopyNumbers(@NotNull final List<FittedRegion> fittedRegions,
            @NotNull final List<PurpleCopyNumber> smoothRegions) {
        return insertSmoothRegions(smoothRegions, fittedRegions);
    }

    @NotNull
    private static List<FittedRegion> insertSmoothRegions(@NotNull final List<PurpleCopyNumber> smoothRegions,
            @NotNull final List<FittedRegion> fittedRegions) {
        BiFunction<PurpleCopyNumber, FittedRegion, FittedRegion> transform =
                (consolidatedRegion, copyNumber) -> ImmutableFittedRegion.builder()
                        .from(copyNumber)
                        .fittedBAF(consolidatedRegion.averageActualBAF())
                        .fittedTumorCopyNumber(consolidatedRegion.averageTumorCopyNumber())
                        .build();

        PurpleRegionZipper zipper = new PurpleRegionZipper(transform);
        RegionZipper.zip(smoothRegions, fittedRegions, zipper);
        return zipper.result();
    }
}
