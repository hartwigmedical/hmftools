package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ExtendNonDiploid extends ExtendRegion
{
    @NotNull
    static List<CombinedRegion> nonDiploid(final List<CombinedRegion> regions)
    {
        return new ExtendNonDiploid().extend(regions);
    }

    private ExtendNonDiploid()
    {
        super(CopyNumberMethod.NON_DIPLOID);
    }

    @NotNull
    @Override
    public List<CombinedRegion> extend(@NotNull final List<CombinedRegion> regions)
    {
        for(int i = 0; i < regions.size(); i++)
        {
            CombinedRegion region = regions.get(i);
            FittedRegion next = i < regions.size() - 1 ? regions.get(i + 1).region() : null;
            if(region.copyNumberMethod().equals(CopyNumberMethod.UNKNOWN) && isEligible(region.region(), next))
            {
                region.setTumorCopyNumber(CopyNumberMethod.NON_DIPLOID, region.region().refNormalisedCopyNumber());
            }
        }
        return super.extend(regions);
    }

    static boolean isEligible(@NotNull final FittedRegion region, @Nullable final FittedRegion neighbour)
    {
        return Doubles.greaterThan(region.observedTumorRatio(), region.observedNormalRatio()) && !region.status()
                .equals(GermlineStatus.DIPLOID) && !region.status().equals(GermlineStatus.NOISE) && region.depthWindowCount() > 0
                && !isBoundByCentromere(region, neighbour) && isBoundBySV(region, neighbour);
    }

    private static boolean isBoundByCentromere(@NotNull final FittedRegion region, @Nullable final FittedRegion neighbour)
    {
        return region.support().equals(SegmentSupport.CENTROMERE) || (neighbour != null && neighbour.support()
                .equals(SegmentSupport.CENTROMERE));
    }

    private static boolean isBoundBySV(@NotNull final FittedRegion region, @Nullable final FittedRegion neighbour)
    {
        return region.support().isSV() || (neighbour != null && neighbour.support().isSV());
    }

    @Override
    protected void extend(final CombinedRegion target, final CombinedRegion neighbour)
    {
        target.extendWithWeightedAverage(neighbour.region());
    }

}
