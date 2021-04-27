package com.hartwig.hmftools.purple.region;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

class GCAccumulator implements Consumer<GCProfile>
{
    private final GenomeRegion mRegion;
    private double mTotalContent;
    private int mCount;

    GCAccumulator(final GenomeRegion region)
    {
        this.mRegion = region;
    }

    double averageGCContent()
    {
        return mCount == 0 ? 0 : mTotalContent / mCount;
    }

    @Override
    public void accept(final GCProfile gcProfile)
    {
        if(gcProfile.isMappable() && gcProfile.start() >= mRegion.start() && gcProfile.end() <= mRegion.end())
        {
            mCount++;
            mTotalContent += gcProfile.gcContent();
        }
    }
}
