package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.Nullable;

class ExtendNonDiploid extends ExtendRegion
{
    static List<CombinedRegion> nonDiploid(final List<CombinedRegion> regions)
    {
        return new ExtendNonDiploid().extend(regions);
    }

    private ExtendNonDiploid()
    {
        super(CopyNumberMethod.NON_DIPLOID);
    }

    @Override
    public List<CombinedRegion> extend(final List<CombinedRegion> regions)
    {
        for(int i = 0; i < regions.size(); i++)
        {
            CombinedRegion region = regions.get(i);
            ObservedRegion next = i < regions.size() - 1 ? regions.get(i + 1).region() : null;
            if(region.copyNumberMethod().equals(CopyNumberMethod.UNKNOWN) && isEligible(region.region(), next))
            {
                region.setTumorCopyNumber(CopyNumberMethod.NON_DIPLOID, region.region().refNormalisedCopyNumber());
            }
        }
        return super.extend(regions);
    }

    static boolean isEligible(final ObservedRegion region, @Nullable final ObservedRegion neighbour)
    {
        return Doubles.greaterThan(region.observedTumorRatio(), region.observedNormalRatio())
                && !region.germlineStatus().equals(GermlineStatus.EXCLUDED)
                && !region.germlineStatus().equals(GermlineStatus.DIPLOID)
                && !region.germlineStatus().equals(GermlineStatus.NOISE)
                && region.depthWindowCount() > 0
                && !isBoundByCentromere(region, neighbour)
                && isBoundBySV(region, neighbour);
    }

    private static boolean isBoundByCentromere(final ObservedRegion region, @Nullable final ObservedRegion neighbour)
    {
        return region.support().equals(SegmentSupport.CENTROMERE) || (neighbour != null && neighbour.support()
                .equals(SegmentSupport.CENTROMERE));
    }

    private static boolean isBoundBySV(final ObservedRegion region, @Nullable final ObservedRegion neighbour)
    {
        return region.support().isSV() || (neighbour != null && neighbour.support().isSV());
    }

    @Override
    protected void extend(final CombinedRegion target, final CombinedRegion neighbour)
    {
        target.extendWithWeightedAverage(neighbour.region());
    }

}
