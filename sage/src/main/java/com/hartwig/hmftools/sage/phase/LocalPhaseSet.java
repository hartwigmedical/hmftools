package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.VariantPhaser.PHASE_BUFFER;

import java.util.Collection;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.SageVariant;

public class LocalPhaseSet extends BufferedPostProcessor
{
    private final PhaseSetCounter mPhaseSetCounter;
    private final Set<Integer> mPassingPhaseSets = Sets.newHashSet();

    public LocalPhaseSet(final PhaseSetCounter phaseSetCounter, final Consumer<SageVariant> consumer)
    {
        super(PHASE_BUFFER, consumer);
        mPhaseSetCounter = phaseSetCounter;
    }

    public Set<Integer> passingPhaseSets()
    {
        return mPassingPhaseSets;
    }

    @Override
    protected void preFlush(final Collection<SageVariant> variants)
    {
        for(SageVariant variant : variants)
        {
            if(variant.isPassing() && variant.hasLocalPhaseSet())
            {
                mPassingPhaseSets.add(variant.localPhaseSet());
            }
        }
    }

    @Override
    protected void processSageVariant(final SageVariant variant, final Collection<SageVariant> variants)
    {
        // now determined earlier

        /*
        final ReadContext newReadContext = variant.readContext();

        for(final SageVariant other : variants)
        {
            final ReadContext otherReadContext = other.readContext();

            if(!rightInLeftDel(other.variant(), variant.variant()))
            {
                int offset = adjustedOffset(other.variant(), variant.variant());
                if(otherReadContext.phased(offset, newReadContext))
                {
                    if(other.localPhaseSet() != 0)
                    {
                        variant.localPhaseSet(other.localPhaseSet());
                    }
                    else if(variant.localPhaseSet() != 0)
                    {
                        other.localPhaseSet(variant.localPhaseSet());
                    }
                    else
                    {
                        int nextLps = mPhaseSetCounter.getNext();
                        other.localPhaseSet(nextLps);
                        variant.localPhaseSet(nextLps);
                    }
                }
            }
        }
        */
    }

    static int positionOffset(final VariantHotspot left, final VariantHotspot right)
    {
        int positionOffset = left.position() - right.position();
        return positionOffset;
    }

    static int adjustedOffset(final VariantHotspot left, final VariantHotspot right)
    {
        int positionOffset = positionOffset(left, right);

        if(positionOffset == 0)
            return 0;

        return positionOffset
                + Math.max(0, left.ref().length() - left.alt().length())
                - Math.max(0, left.alt().length() - left.ref().length());
    }

    static boolean rightInLeftDel(final VariantHotspot left, final VariantHotspot right)
    {
        if(left.ref().length() > left.alt().length())
        {
            int deleteEnd = left.position() + left.ref().length() - 1;
            return right.position() > left.position() && right.position() <= deleteEnd;
        }

        return false;
    }
}
