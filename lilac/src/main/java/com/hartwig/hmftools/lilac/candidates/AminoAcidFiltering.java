package com.hartwig.hmftools.lilac.candidates;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilackt.SequenceCount;
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci;

import java.util.List;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class AminoAcidFiltering
{
    private final Set<Integer> mAminoAcidBoundaries;

    public AminoAcidFiltering(@NotNull Set<Integer> aminoAcidBoundaries)
    {
        mAminoAcidBoundaries = aminoAcidBoundaries;
    }

    public final List<HlaSequenceLoci> aminoAcidCandidates(@NotNull List<HlaSequenceLoci> candidates, @NotNull SequenceCount aminoAcidCount)
    {
        // TODO
        return Lists.newArrayList();

        /*
        List result = candidates;
        Set locations =
                CollectionsKt.subtract((Iterable) CollectionsKt.toSet((Iterable) ((Iterable) RangesKt.until((int) 0, (int) aminoAcidCount.getLength()))), (Iterable) this.mAminoAcidBoundaries);
        Iterator iterator = locations.iterator();
        while(iterator.hasNext())
        {
            void $receiver$iv$iv;
            Iterable $receiver$iv;
            int loci = ((Number) iterator.next()).intValue();
            Collection<String> expectedSequences = aminoAcidCount.sequenceAt(loci);
            Iterable iterable = $receiver$iv = (Iterable) result;
            Collection destination$iv$iv = new ArrayList();
            for(Object element$iv$iv : $receiver$iv$iv)
            {
                HlaSequenceLoci it = (HlaSequenceLoci) element$iv$iv;
                boolean bl = false;
                if(!it.consistentWithAny(expectedSequences, loci))
                {
                    continue;
                }
                destination$iv$iv.add(element$iv$iv);
            }
            result = (List) destination$iv$iv;
        }
        return result;

         */
    }

}
