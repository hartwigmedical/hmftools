package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.Phase.PHASE_BUFFER;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

class LocalRealignSet extends BufferedPostProcessor
{

    private int realign;

    LocalRealignSet(@NotNull final Consumer<SageVariant> consumer)
    {
        super(PHASE_BUFFER, consumer);
    }

    @Override
    protected void processSageVariant(@NotNull final SageVariant newEntry, @NotNull final Collection<SageVariant> buffer)
    {
        final ReadContext newReadContext = newEntry.readContext();
        boolean newIsIndel = newEntry.isIndel();

        for(final SageVariant oldEntry : buffer)
        {
            boolean oldIsIndel = oldEntry.isIndel();
            if(newIsIndel || oldIsIndel)
            {
                final ReadContext oldReadContext = oldEntry.readContext();
                int positionOffset = LocalPhaseSet.positionOffset(oldEntry.variant(), newEntry.variant());
                int offset = LocalPhaseSet.adjustedOffset(oldEntry.variant(), newEntry.variant());

                if(positionOffset != offset && oldReadContext.phased(positionOffset, newReadContext))
                {
                    if(oldEntry.localRealignSet() != 0)
                    {
                        newEntry.localRealignSet(oldEntry.localRealignSet());
                    }
                    else if(newEntry.localRealignSet() != 0)
                    {
                        oldEntry.localRealignSet(newEntry.localRealignSet());
                    }
                    else
                    {
                        realign++;
                        oldEntry.localRealignSet(realign);
                        newEntry.localRealignSet(realign);
                    }
                }
            }
        }
    }
}
