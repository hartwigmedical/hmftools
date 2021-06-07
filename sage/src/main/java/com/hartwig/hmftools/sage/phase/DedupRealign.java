package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.Phase.PHASE_BUFFER;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.jetbrains.annotations.NotNull;

public class DedupRealign extends BufferedPostProcessor
{

    public DedupRealign(final Consumer<SageVariant> consumer)
    {
        super(PHASE_BUFFER, consumer);
    }

    @Override
    protected void processSageVariant(@NotNull final SageVariant newVariant, @NotNull final Collection<SageVariant> buffer)
    {
        // Empty
    }

    @Override
    protected void preFlush(@NotNull final Collection<SageVariant> variants)
    {
        super.preFlush(variants);
        final Set<Integer> localRealignedSets = variants.stream()
                .filter(x -> x.filters().isEmpty())
                .map(SageVariant::localRealignSet)
                .filter(x -> x > 0)
                .collect(Collectors.toSet());

        for(Integer localRealignedSet : localRealignedSets)
        {
            process(localRealignedSet, variants);
        }
    }

    private void process(int localRealignSet, @NotNull final Collection<SageVariant> buffer)
    {
        final List<SageVariant> localRealigned = buffer.stream()
                .filter(x -> x.filters().isEmpty())
                .filter(x -> x.localRealignSet() == localRealignSet)
                .collect(Collectors.toList());

        if(localRealigned.size() > 1)
        {
            SageVariant keeper = localRealigned.get(0);
            for(SageVariant sageVariant : localRealigned)
            {
                if(sageVariant.totalQuality() > keeper.totalQuality())
                {
                    keeper = sageVariant;
                }
            }

            for(SageVariant sageVariant : localRealigned)
            {
                boolean keep =
                        sageVariant.equals(keeper) || (keeper.localPhaseSet() > 0 && keeper.localPhaseSet() == sageVariant.localPhaseSet());
                if(!keep)
                {
                    sageVariant.filters().add(SageVCF.DEDUP_FILTER);
                }
            }
        }
    }
}
