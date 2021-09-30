package com.hartwig.hmftools.purple.gene;

import java.util.List;
import java.util.function.BiFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PurpleRegionZipper implements RegionZipperHandler<PurpleCopyNumber, FittedRegion>
{
    private final BiFunction<PurpleCopyNumber, FittedRegion, FittedRegion> mTransform;
    private final List<FittedRegion> mResult = Lists.newArrayList();

    @Nullable
    private PurpleCopyNumber mConsolidatedRegion;

    private PurpleRegionZipper(@NotNull final BiFunction<PurpleCopyNumber, FittedRegion, FittedRegion> transform)
    {
        mTransform = transform;
    }

    @NotNull
    private List<FittedRegion> result()
    {
        return mResult;
    }

    @Override
    public void enterChromosome(@NotNull final String chromosome)
    {
        // Empty
    }

    @Override
    public void primary(@NotNull final PurpleCopyNumber region)
    {
        mConsolidatedRegion = region;
    }

    @Override
    public void secondary(@NotNull final FittedRegion copyNumber)
    {
        if(mConsolidatedRegion != null && mConsolidatedRegion.overlaps(copyNumber))
        {
            mResult.add(mTransform.apply(mConsolidatedRegion, copyNumber));
        }
        else
        {
            mResult.add(copyNumber);
        }
    }

    public static List<FittedRegion> updateRegionsWithCopyNumbers(@NotNull final List<FittedRegion> fittedRegions,
            @NotNull final List<PurpleCopyNumber> smoothRegions)
    {
        return insertSmoothRegions(smoothRegions, fittedRegions);
    }

    @NotNull
    private static List<FittedRegion> insertSmoothRegions(@NotNull final List<PurpleCopyNumber> smoothRegions,
            @NotNull final List<FittedRegion> fittedRegions)
    {
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
