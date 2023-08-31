package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public final class ExtendUtils
{
    static boolean doNotExtend(final CombinedRegion target, final ObservedRegion neighbour)
    {
        return breakForCentromereStart(target, neighbour) || breakForStructuralVariant(target, neighbour);
    }

    private static boolean breakForCentromereStart(final CombinedRegion target, final ObservedRegion neighbour)
    {
        if(target.start() < neighbour.start())
        {
            return neighbour.support() == SegmentSupport.CENTROMERE;
        }

        return target.region().support() == SegmentSupport.CENTROMERE;
    }

    private static boolean breakForStructuralVariant(final CombinedRegion target, final ObservedRegion neighbour)
    {
        if(target.start() < neighbour.start())
        {
            return neighbour.support().isSV();
        }

        return target.region().support().isSV();
    }

    public static List<CombinedRegion> populateUnknown(final List<CombinedRegion> regions, final CobaltChromosomes cobaltChromosomes)
    {
        for(int i = 0; i < regions.size(); i++)
        {
            final CombinedRegion region = regions.get(i);

            if(region.copyNumberMethod() == CopyNumberMethod.UNKNOWN)
            {
                double normalCopyNumber = 2 * cobaltChromosomes.get(region.chromosome()).actualRatio();
                region.setTumorCopyNumber(CopyNumberMethod.UNKNOWN, normalCopyNumber);

                if(region.support() == SegmentSupport.NONE && i > 0)
                {
                    final CombinedRegion prev = regions.get(i - 1);
                    if(prev.copyNumberMethod() == CopyNumberMethod.UNKNOWN)
                    {
                        prev.extend(region.region());
                        regions.remove(i);
                        i--;
                    }
                }
            }
        }

        return regions;
    }

}
