package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.jetbrains.annotations.NotNull;

public class DedupMnv extends BufferedPostProcessor
{
    private static final int BUFFER = 10;

    public DedupMnv(final Consumer<SageVariant> consumer)
    {
        super(BUFFER, consumer);
    }

    @Override
    protected void processSageVariant(@NotNull final SageVariant newVariant, @NotNull final Collection<SageVariant> buffer)
    {
        int lps = newVariant.localPhaseSet();

        if(newVariant.isPassing() && !newVariant.isIndel() && lps > 0)
        {

            int newVariantSize = newVariant.alt().length();
            for(final SageVariant oldVariant : buffer)
            {
                if(oldVariant.isPassing() && !oldVariant.isIndel() && oldVariant.localPhaseSet() == lps)
                {
                    int oldVariantSize = oldVariant.alt().length();
                    if(newVariantSize != oldVariantSize)
                    {
                        final SageVariant shorter;
                        final SageVariant longer;
                        if(newVariantSize > oldVariantSize)
                        {
                            shorter = oldVariant;
                            longer = newVariant;
                        }
                        else
                        {
                            shorter = newVariant;
                            longer = oldVariant;
                        }
                        if(longerContainsShorter(shorter, longer))
                        {
                            shorter.filters().add(SageVCF.DEDUP_FILTER);
                        }
                    }
                }
            }
        }
    }
}
