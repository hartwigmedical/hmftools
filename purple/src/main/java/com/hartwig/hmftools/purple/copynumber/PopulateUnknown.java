package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public class PopulateUnknown
{
    private final CobaltChromosomes mCobaltChromosomes;

    public PopulateUnknown(final CobaltChromosomes cobaltChromosomes)
    {
        mCobaltChromosomes = cobaltChromosomes;
    }

    public List<CombinedRegion> populateUnknown(@NotNull final List<CombinedRegion> regions)
    {
        for(int i = 0; i < regions.size(); i++)
        {
            final CombinedRegion region = regions.get(i);

            if(region.copyNumberMethod() == CopyNumberMethod.UNKNOWN)
            {
                double normalCopyNumber = 2 * mCobaltChromosomes.get(region.chromosome()).actualRatio();
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