package com.hartwig.hmftools.sage.phase;

import java.util.Collection;
import java.util.function.Consumer;

import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import org.jetbrains.annotations.NotNull;

public class DedupIndel extends BufferedPostProcessor
{
    DedupIndel(final Consumer<SageVariant> consumer)
    {
        super(0, consumer);
    }

    @Override
    protected void processSageVariant(@NotNull final SageVariant variant, @NotNull final Collection<SageVariant> buffer)
    {
        if(isPassingPhasedIndel(variant))
        {
            for(final SageVariant other : buffer)
            {
                if(isPassingPhasedIndel(other) && variant.localPhaseSet() == other.localPhaseSet())
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
    }

    private void processDel(@NotNull final SageVariant left, @NotNull final SageVariant right)
    {
        if(!left.alt().equals(right.alt()))
        {
            return;
        }

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

    private void processIns(@NotNull final SageVariant left, @NotNull final SageVariant right)
    {
        if(!left.ref().equals(right.ref()))
        {
            return;
        }

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

    private static boolean isPassingPhasedIndel(@NotNull final SageVariant newEntry)
    {
        return newEntry.isPassing() && newEntry.localPhaseSet() > 0 && newEntry.isIndel();
    }
}
