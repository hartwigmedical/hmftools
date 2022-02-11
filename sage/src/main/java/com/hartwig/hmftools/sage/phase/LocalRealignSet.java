package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.VariantDeduperOld.PHASE_BUFFER;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
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
            if(!other.isIndel() && !newIsIndel)
                continue;

            if(!variant.hasLocalPhaseSets() || !other.hasMatchingLps(variant.localPhaseSets()))
                continue;

            int positionOffset = positionOffset(other.variant(), variant.variant());
            int offset = adjustedOffset(other.variant(), variant.variant());

            if(positionOffset != offset && other.readContext().phased(positionOffset, newReadContext))
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

    public static int positionOffset(final VariantHotspot left, final VariantHotspot right)
    {
        int positionOffset = left.position() - right.position();
        return positionOffset;
    }

    public static int adjustedOffset(final VariantHotspot left, final VariantHotspot right)
    {
        int positionOffset = positionOffset(left, right);

        if(positionOffset == 0)
            return 0;

        return positionOffset
                + Math.max(0, left.ref().length() - left.alt().length())
                - Math.max(0, left.alt().length() - left.ref().length());
    }

}
