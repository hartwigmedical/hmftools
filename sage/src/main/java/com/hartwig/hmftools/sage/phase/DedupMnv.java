package com.hartwig.hmftools.sage.phase;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

public class DedupMnv extends BufferedPostProcessor
{
    private static final int BUFFER = 10;

    public DedupMnv(final Consumer<SageVariant> consumer)
    {
        super(BUFFER, consumer);
    }

    @Override
    protected void processSageVariant(final SageVariant variant, final Collection<SageVariant> variants)
    {
        if(!variant.isPassing() || variant.isIndel() || variant.localPhaseSet() <= 0)
            return;

        int lps = variant.localPhaseSet();
        int newVariantSize = variant.alt().length();

        for(final SageVariant other : variants)
        {
            if(other.isPassing() && !other.isIndel() && other.localPhaseSet() == lps)
            {
                int oldVariantSize = other.alt().length();
                if(newVariantSize != oldVariantSize)
                {
                    final SageVariant shorter;
                    final SageVariant longer;
                    if(newVariantSize > oldVariantSize)
                    {
                        shorter = other;
                        longer = variant;
                    }
                    else
                    {
                        shorter = variant;
                        longer = other;
                    }
                    if(longerContainsShorter(shorter, longer))
                    {
                        shorter.filters().add(VariantVCF.DEDUP_FILTER);
                    }
                }
            }
        }
    }
}
