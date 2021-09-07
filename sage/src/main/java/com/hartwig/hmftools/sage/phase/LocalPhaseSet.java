package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.Phase.PHASE_BUFFER;

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
    private int mPhase;
    private final Set<Integer> mPassingPhaseSets = Sets.newHashSet();

    public LocalPhaseSet(@NotNull final Consumer<SageVariant> consumer)
    {
        super(PHASE_BUFFER, consumer);
    }

    @NotNull
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
    protected void processSageVariant(final SageVariant newEntry, final Collection<SageVariant> buffer)
    {
        final ReadContext newReadContext = newEntry.readContext();
        for(final SageVariant oldEntry : buffer)
        {
            final ReadContext oldReadContext = oldEntry.readContext();

            if(!rightInLeftDel(oldEntry.variant(), newEntry.variant()))
            {
                int offset = adjustedOffset(oldEntry.variant(), newEntry.variant());
                if(oldReadContext.phased(offset, newReadContext))
                {
                    if(oldEntry.localPhaseSet() != 0)
                    {
                        newEntry.localPhaseSet(oldEntry.localPhaseSet());
                    }
                    else if(newEntry.localPhaseSet() != 0)
                    {
                        oldEntry.localPhaseSet(newEntry.localPhaseSet());
                    }
                    else
                    {
                        mPhase++;
                        oldEntry.localPhaseSet(mPhase);
                        newEntry.localPhaseSet(mPhase);
                    }
                }
            }
        }
    }

    static int positionOffset(final VariantHotspot left, final VariantHotspot right)
    {
        long positionOffset = left.position() - right.position();
        return (int) (positionOffset);
    }

    static int adjustedOffset(final VariantHotspot left, final VariantHotspot right)
    {
        long positionOffset = positionOffset(left, right);
        if(positionOffset == 0)
        {
            return 0;
        }

        return (int) (positionOffset + Math.max(0, left.ref().length() - left.alt().length()) - Math.max(0,
                left.alt().length() - left.ref().length()));
    }

    static boolean rightInLeftDel(final VariantHotspot left, final VariantHotspot right)
    {
        if(left.ref().length() > left.alt().length())
        {
            long deleteEnd = left.position() + left.ref().length() - 1;
            return right.position() > left.position() && right.position() <= deleteEnd;
        }

        return false;
    }
}
