package com.hartwig.hmftools.lilac.candidates;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragment;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import java.util.List;
import java.util.Set;

public class NucleotideFiltering
{
    private final int mMinNucleotideCount;
    private final Set<Integer> mAminoAcidBoundaries;

    public NucleotideFiltering(int minNucleotideCount, final Set<Integer> aminoAcidBoundaries)
    {
        mMinNucleotideCount = minNucleotideCount;
        mAminoAcidBoundaries = aminoAcidBoundaries;
    }

    public final List<HlaSequenceLoci> filterCandidatesOnAminoAcidBoundaries(
            final List<HlaSequenceLoci> candidates, final List<NucleotideFragment> fragments)
    {
        final List<HlaSequenceLoci> results = Lists.newArrayList();

        for(int boundary : mAminoAcidBoundaries)
        {
            int nucleotideStart = boundary * 3;
            final List<String> startSequences = nucleotideSequence(fragments, Sets.newHashSet(nucleotideStart));
            final List<String> endSequences = nucleotideSequence(fragments, Sets.newHashSet(nucleotideStart + 1, nucleotideStart + 2));

            // CHECK
            // result = result.filter { it.consistentWithAny(nucleotideStart, startSequences, endSequences) }
        }

        return results;
    }

    private final boolean consistentWithAny(
            final HlaSequenceLoci seqLoci, int startLoci, final List<String> startSequences, final List<String> endSequences)
    {
        return seqLoci.consistentWithAny(startSequences, Sets.newHashSet(startLoci))
            && seqLoci.consistentWithAny(endSequences, Sets.newHashSet(startLoci + 1, startLoci + 2));
    }

    private final List<String> nucleotideSequence(List<NucleotideFragment> fragments, final Set<Integer> nucleotideIndices)
    {
        return Lists.newArrayList();

        /*
        String string;
        Object object;
        NucleotideFragment it;
        Object $receiver$iv$iv;
        Object $receiver$iv = fragments;
        Iterable iterable = $receiver$iv;
        Object destination$iv$iv = new ArrayList();
        Object object2 = $receiver$iv$iv.iterator();
        while(object2.hasNext())
        {
            Object element$iv$iv = object2.next();
            it = (NucleotideFragment) element$iv$iv;
            boolean bl = false;
            if(!it.containsAllNucleotides(Arrays.copyOf(nucleotideIndices, nucleotideIndices.length)))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        object2 = $receiver$iv$iv.iterator();
        while(object2.hasNext())
        {
            Object item$iv$iv = object2.next();
            it = (NucleotideFragment) item$iv$iv;
            object = destination$iv$iv;
            boolean bl = false;
            string = it.nucleotides(Arrays.copyOf(nucleotideIndices, nucleotideIndices.length));
            object.add(string);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv = GroupingKt.eachCount((Grouping) ((Grouping) new Grouping<String, String>((Iterable) $receiver$iv)
        {
            final Iterable receiver$0;

            {
                receiver$0 = $receiver;
            }

            @NotNull
            public Iterator<String> sourceIterator()
            {
                return receiver$0.iterator();
            }

            public Object keyOf(Object element)
            {
                void var2_2;
                String it = (String) element;
                boolean bl = false;
                return var2_2;
            }
        }));
        Object $i$f$groupingBy = $receiver$iv;
        destination$iv$iv = new LinkedHashMap();
        object2 = $receiver$iv$iv;
        Iterator iterator = object2.entrySet().iterator();
        while(iterator.hasNext())
        {
            Map.Entry element$iv$iv;
            Map.Entry it2 = element$iv$iv = iterator.next();
            boolean bl = false;
            if(!(((Number) it2.getValue()).intValue() >= mMinNucleotideCount))
            {
                continue;
            }
            destination$iv$iv.put(element$iv$iv.getKey(), element$iv$iv.getValue());
        }
        $receiver$iv = destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList($receiver$iv.size());
        object2 = $receiver$iv$iv;
        iterator = object2.entrySet().iterator();
        while(iterator.hasNext())
        {
            Map.Entry item$iv$iv;
            Map.Entry it2 = item$iv$iv = iterator.next();
            object = destination$iv$iv;
            boolean bl = false;
            string = (String) it2.getKey();
            object.add(string);
        }
        return (List) destination$iv$iv;
        */
    }

}
