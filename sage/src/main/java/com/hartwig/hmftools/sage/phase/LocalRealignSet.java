package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.VariantPhaser.PHASE_BUFFER;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.SageVariant;

class LocalRealignSet extends BufferedPostProcessor
{
    private int mRealign;

    LocalRealignSet(final Consumer<SageVariant> consumer)
    {
        super(PHASE_BUFFER, consumer);
    }

    @Override
    protected void processSageVariant(final SageVariant variant, final Collection<SageVariant> variants)
    {
        final ReadContext newReadContext = variant.readContext();
        boolean newIsIndel = variant.isIndel();

        for(final SageVariant other : variants)
        {
            boolean oldIsIndel = other.isIndel();
            if(newIsIndel || oldIsIndel)
            {
                final ReadContext oldReadContext = other.readContext();
                int positionOffset = LocalPhaseSet.positionOffset(other.variant(), variant.variant());
                int offset = LocalPhaseSet.adjustedOffset(other.variant(), variant.variant());

                if(positionOffset != offset && oldReadContext.phased(positionOffset, newReadContext))
                {
                    if(other.localRealignSet() != 0)
                    {
                        variant.localRealignSet(other.localRealignSet());
                    }
                    else if(variant.localRealignSet() != 0)
                    {
                        other.localRealignSet(variant.localRealignSet());
                    }
                    else
                    {
                        mRealign++;
                        other.localRealignSet(mRealign);
                        variant.localRealignSet(mRealign);
                    }
                }
            }
        }
    }
}
