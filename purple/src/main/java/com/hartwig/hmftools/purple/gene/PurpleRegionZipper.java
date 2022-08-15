package com.hartwig.hmftools.purple.gene;

import java.util.List;
import java.util.function.BiFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PurpleRegionZipper implements RegionZipperHandler<PurpleCopyNumber, ObservedRegion>
{
    private final BiFunction<PurpleCopyNumber, ObservedRegion, ObservedRegion> mTransform;
    private final List<ObservedRegion> mResult = Lists.newArrayList();

    @Nullable
    private PurpleCopyNumber mConsolidatedRegion;

    private PurpleRegionZipper(final BiFunction<PurpleCopyNumber, ObservedRegion, ObservedRegion> transform)
    {
        mTransform = transform;
    }

    @NotNull
    private List<ObservedRegion> result()
    {
        return mResult;
    }

    @Override
    public void enterChromosome(final String chromosome)
    {
        // Empty
    }

    @Override
    public void primary(final PurpleCopyNumber region)
    {
        mConsolidatedRegion = region;
    }

    @Override
    public void secondary(final ObservedRegion copyNumber)
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

    public static List<ObservedRegion> updateRegionsWithCopyNumbers(
            final List<ObservedRegion> fittedRegions, final List<PurpleCopyNumber> smoothRegions)
    {
        return insertSmoothRegions(smoothRegions, fittedRegions);
    }

    @NotNull
    private static List<ObservedRegion> insertSmoothRegions(final List<PurpleCopyNumber> smoothRegions,
            final List<ObservedRegion> fittedRegions)
    {
        BiFunction<PurpleCopyNumber, ObservedRegion, ObservedRegion> transform =
                (consolidatedRegion, copyNumber) -> transformRegion(
                        copyNumber, consolidatedRegion.averageActualBAF(), consolidatedRegion.averageTumorCopyNumber());

        PurpleRegionZipper zipper = new PurpleRegionZipper(transform);
        RegionZipper.zip(smoothRegions, fittedRegions, zipper);
        return zipper.result();
    }

    private static ObservedRegion transformRegion(final ObservedRegion region, double fittedBAF, double fittedTumorCopyNumber)
    {
        ObservedRegion newRegion = ObservedRegion.from(region);
        newRegion.setFittedBAF(fittedBAF);
        newRegion.setFittedTumorCopyNumber(fittedTumorCopyNumber);
        return newRegion;
    }
}
