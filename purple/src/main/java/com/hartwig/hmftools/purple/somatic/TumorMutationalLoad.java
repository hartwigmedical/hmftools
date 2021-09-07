package com.hartwig.hmftools.purple.somatic;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class TumorMutationalLoad
{
    private static final VariantContextFilter PASS = new PassingVariantFilter();

    private int mLoad;
    private int mBurden;

    public int load() {
        return mLoad;
    }

    public double burdenPerMb()
    {
        return mBurden / MicrosatelliteIndels.NUMBER_OF_MB_PER_GENOME;
    }

    public void processVariant(VariantContextDecorator variant)
    {
        if(variant.isPass())
        {
            mBurden++;

            final VariantImpact variantImpact = variant.variantImpact();

            if (variantImpact.WorstCodingEffect.equals(CodingEffect.MISSENSE))
                mLoad++;
        }
    }
}
