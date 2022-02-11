package com.hartwig.hmftools.sage.phase;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

public class DedupIndelOld extends BufferedPostProcessor
{
    DedupIndelOld(final Consumer<SageVariant> consumer)
    {
        super(0, consumer);
    }

    @Override
    protected void processSageVariant(final SageVariant variant, final Collection<SageVariant> variants)
    {
        if(!isPassingPhasedIndel(variant))
            return;

        for(final SageVariant other : variants)
        {
            if(isPassingPhasedIndel(other) && other.hasMatchingLps(variant.localPhaseSets()))
            {
                if(variant.isDelete() && other.isDelete())
                {
                    processDel(variant, other);
                }

                if(variant.isInsert() && other.isInsert())
                {
                    processIns(variant, other);
                }
            }
        }
    }

    private void processDel(final SageVariant left, final SageVariant right)
    {
        if(!left.alt().equals(right.alt()))
            return;

        final SageVariant shorter;
        final SageVariant longer;
        if(left.ref().length() < right.ref().length())
        {
            shorter = left;
            longer = right;
        }
        else
        {
            shorter = right;
            longer = left;
        }

        if(longer.ref().substring(0, shorter.ref().length()).equals(shorter.ref()))
        {
            longer.filters().add(VariantVCF.DEDUP_FILTER);
        }
    }

    private void processIns(final SageVariant left, final SageVariant right)
    {
        if(!left.ref().equals(right.ref()))
            return;

        final SageVariant shorter;
        final SageVariant longer;
        if(left.alt().length() < right.alt().length())
        {
            shorter = left;
            longer = right;
        }
        else
        {
            shorter = right;
            longer = left;
        }

        if(longer.alt().substring(0, shorter.alt().length()).equals(shorter.alt()))
        {
            longer.filters().add(VariantVCF.DEDUP_FILTER);
        }
    }

    private static boolean isPassingPhasedIndel(final SageVariant variant)
    {
        return variant.isPassing() && variant.hasLocalPhaseSets() && variant.isIndel();
    }
}
