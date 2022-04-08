package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.purple.config.PurpleConstants.MB_PER_GENOME;

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

    public int load()
    {
        return mLoad;
    }
    public int burden() { return mBurden; }

    public int tml() { return mTargetRegions.calcTml(mLoad); }

    public double burdenPerMb() { return mBurden / MB_PER_GENOME; }

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
