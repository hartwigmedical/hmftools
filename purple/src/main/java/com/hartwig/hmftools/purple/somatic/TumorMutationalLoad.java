package com.hartwig.hmftools.purple.somatic;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.purple.config.TargetRegionsData;

public class TumorMutationalLoad
{
    private final TargetRegionsData mTargetRegions;
    private int mLoad;
    private int mBurden;

    public TumorMutationalLoad(final TargetRegionsData targetRegions)
    {
        mTargetRegions = targetRegions;
        mLoad = 0;
        mBurden = 0;
    }

    public int load() { return mTargetRegions.calcTml(mLoad); }

    public double burdenPerMb()
    {
        return mTargetRegions.calcTmb(mBurden);
    }

    public void processVariant(final SomaticVariant variant)
    {
        if(mTargetRegions.hasTargetRegions() && !mTargetRegions.inTargetRegions(variant.chromosome(), variant.position()))
            return;

        mBurden++;

        final VariantImpact variantImpact = variant.variantImpact();

        if(variantImpact.WorstCodingEffect.equals(CodingEffect.MISSENSE))
            mLoad++;
    }
}
