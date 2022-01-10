package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.VariantPhaser.PHASE_BUFFER;

import java.util.Collection;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

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
        super.preFlush(variants);
        for(SageVariant variant : variants)
        {
            if(variant.isPassing() && variant.localPhaseSet() > 0)
            {
                mPassingPhaseSets.add(variant.localPhaseSet());
            }
        }
    }

    @Override
    protected void processSageVariant(final SageVariant newVariant, final Collection<SageVariant> existingVariants)
    {
        final ReadContext newReadContext = newVariant.readContext();

        for(final SageVariant existingVariant : existingVariants)
        {
            final ReadContext oldReadContext = existingVariant.readContext();

            if(!rightInLeftDel(existingVariant.variant(), newVariant.variant()))
            {
                int offset = adjustedOffset(existingVariant.variant(), newVariant.variant());
                if(oldReadContext.phased(offset, newReadContext))
                {
                    if(existingVariant.localPhaseSet() != 0)
                    {
                        newVariant.localPhaseSet(existingVariant.localPhaseSet());
                    }
                    else if(newVariant.localPhaseSet() != 0)
                    {
                        existingVariant.localPhaseSet(newVariant.localPhaseSet());
                    }
                    else
                    {
                        int nextLps = mPhaseSetCounter.getNext();
                        existingVariant.localPhaseSet(nextLps);
                        newVariant.localPhaseSet(nextLps);
                    }
                }
            }
        }
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
