package com.hartwig.hmftools.sage.phase;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class MixedSomaticGermlineIdentifier extends BufferedPostProcessor
{
    private static final int MAX_DISTANCE = 10;

    private int mCombinedImpact;

    public MixedSomaticGermlineIdentifier(final Consumer<SageVariant> consumer)
    {
        super(MAX_DISTANCE, consumer);
    }

    @Override
    protected void processSageVariant(@NotNull final SageVariant newVariant, @NotNull final Collection<SageVariant> buffer)
    {
        if(!newVariant.isIndel())
        {
            boolean newIsPassingMnv = isPassingMnv(newVariant);
            boolean newIsGermlineSnv = isGermlineFilteredSnv(newVariant);
            if(newIsPassingMnv || newIsGermlineSnv)
            {
                for(final SageVariant other : buffer)
                {
                    if(newIsPassingMnv && isGermlineFilteredSnv(other))
                    {
                        process(newVariant, other);
                    }
                    else if(newIsGermlineSnv && isPassingMnv(other))
                    {
                        process(other, newVariant);
                    }
                }
            }
        }
    }

    private void process(@NotNull final SageVariant mnv, @NotNull final SageVariant germlineSnv)
    {
        if(DedupMnv.longerContainsShorter(germlineSnv, mnv))
        {
            if(mnv.mixedGermlineImpact() == 0)
            {
                mnv.mixedGermlineImpact(++mCombinedImpact);
            }
            germlineSnv.mixedGermlineImpact(mnv.mixedGermlineImpact());
        }
    }

    private static boolean isGermlineFilteredSnv(final SageVariant newVariant)
    {
        return !newVariant.isPassing() && newVariant.isSnv() && SoftFilter.isGermlineAndNotTumorFiltered(newVariant.filters());
    }

    private static boolean isPassingMnv(final SageVariant newVariant)
    {
        return newVariant.isPassing() && newVariant.isMnv();
    }

}
