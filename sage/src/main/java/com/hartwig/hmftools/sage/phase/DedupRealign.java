package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.VariantDeduperOld.PHASE_BUFFER;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

public class DedupRealign extends BufferedPostProcessor
{
    public DedupRealign(final Consumer<SageVariant> consumer)
    {
        super(PHASE_BUFFER, consumer);
    }

    @Override
    protected void processSageVariant(final SageVariant variant, final Collection<SageVariant> variants)
    {
        // Empty
    }

    @Override
    protected void preFlush(final Collection<SageVariant> variants)
    {
        final Set<Integer> localRealignedSets = variants.stream()
                .filter(x -> x.filters().isEmpty())
                .filter(x -> x.hasLocalRealignSet())
                .map(SageVariant::localRealignSet)
                .collect(Collectors.toSet());

        for(Integer localRealignedSet : localRealignedSets)
        {
            process(localRealignedSet, variants);
        }
    }

    private void process(int localRealignSet, final Collection<SageVariant> buffer)
    {
        final List<SageVariant> localRealigned = buffer.stream()
                .filter(x -> x.filters().isEmpty())
                .filter(x -> x.localRealignSet() == localRealignSet)
                .collect(Collectors.toList());

        if(localRealigned.size() > 1)
        {
            SageVariant keeper = localRealigned.get(0);
            for(SageVariant variant : localRealigned)
            {
                if(variant.totalQuality() > keeper.totalQuality())
                {
                    keeper = variant;
                }
            }

            for(SageVariant variant : localRealigned)
            {
                boolean keep = variant.equals(keeper)
                        || (keeper.hasLocalPhaseSets() && keeper.hasMatchingLps(variant.localPhaseSets()));

                if(!keep)
                {
                    variant.filters().add(VariantVCF.DEDUP_FILTER);
                }
            }
        }
    }
}
