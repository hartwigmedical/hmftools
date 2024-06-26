package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.sage.filter.SoftFilter.DEDUP_MATCH;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public final class DedupMatching
{
    public static void dedupMatchingVariants(final List<SageVariant> variants)
    {
        // remove one of the variants if they are identical but for read-context
        List<SageVariant> passingVariants = variants.stream()
                .filter(x -> x.isPassing())
                .collect(Collectors.toList());;

        int index = 0;
        while(index < passingVariants.size() - 1)
        {
            SageVariant firstVariant = passingVariants.get(index);

            // look for an exact match
            int nextIndex = index + 1;
            SageVariant nextVariant = passingVariants.get(nextIndex);

            if(exactMatch(firstVariant, nextVariant))
            {
                // keep the variant with higher qual and merge the local phase sets
                if(firstVariant.totalQuality() >= nextVariant.totalQuality())
                {
                    addLocalPhaseSets(nextVariant, firstVariant);
                    nextVariant.filters().add(DEDUP_MATCH);
                }
                else
                {
                    addLocalPhaseSets(firstVariant, nextVariant);
                    firstVariant.filters().add(DEDUP_MATCH);
                }

                index = nextIndex + 1;
            }
            else
            {
                index = nextIndex;
            }
        }
    }

    private static boolean exactMatch(final SageVariant first, final SageVariant second)
    {
        return first.position() == second.position() && first.alt().equals(second.alt()) && first.ref().equals(second.ref());
    }

    private static void addLocalPhaseSets(final SageVariant source, final SageVariant dest)
    {
        if(!source.hasLocalPhaseSets() || dest.tumorReadCounters().isEmpty())
            return;

        List<Integer> localPhaseSets = source.localPhaseSets();
        List<Integer> lpsCounts = source.localPhaseSetCounts();
        ReadContextCounter destReadCounter = dest.tumorReadCounters().get(0);

        for(int i = 0; i < localPhaseSets.size(); ++i)
        {
            destReadCounter.addLocalPhaseSet(localPhaseSets.get(i), lpsCounts.get(i), 0);
        }
    }
}
