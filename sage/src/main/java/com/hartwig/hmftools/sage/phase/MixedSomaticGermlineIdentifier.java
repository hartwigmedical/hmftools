package com.hartwig.hmftools.sage.phase;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.common.SageVariant;

public class MixedSomaticGermlineIdentifier extends BufferedPostProcessor
{
    private static final int MAX_DISTANCE = 10;

    private int mCombinedImpact;

    public MixedSomaticGermlineIdentifier(final Consumer<SageVariant> consumer)
    {
        super(MAX_DISTANCE, consumer);
    }

    @Override
    protected void processSageVariant(final SageVariant variant, final Collection<SageVariant> variants)
    {
        if(variant.isIndel())
            return;

        boolean newIsPassingMnv = isPassingMnv(variant);
        boolean newIsGermlineSnv = isGermlineFilteredSnv(variant);
        if(newIsPassingMnv || newIsGermlineSnv)
        {
            for(final SageVariant other : variants)
            {
                if(newIsPassingMnv && isGermlineFilteredSnv(other))
                {
                    process(variant, other);
                }
                else if(newIsGermlineSnv && isPassingMnv(other))
                {
                    process(other, variant);
                }
            }
        }
    }

    private void process(final SageVariant mnv, final SageVariant germlineSnv)
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
