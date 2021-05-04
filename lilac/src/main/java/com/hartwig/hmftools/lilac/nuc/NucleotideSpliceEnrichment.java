package com.hartwig.hmftools.lilac.nuc;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.SequenceCount;

import java.util.List;
import java.util.Set;

public class NucleotideSpliceEnrichment
{
    private final int mMinBaseQuality;
    private final int mMinBaseCount;
    private final Set<Integer> mAminoAcidBoundary;

    public NucleotideSpliceEnrichment(int minBaseQuality, int minBaseCount, final Set<Integer> aminoAcidBoundary)
    {
        mMinBaseQuality = minBaseQuality;
        mMinBaseCount = minBaseCount;
        mAminoAcidBoundary = aminoAcidBoundary;
    }

    public final List<NucleotideFragment> enrich(final List<NucleotideFragment> fragments)
    {
        // TODO
        return Lists.newArrayList();
        /*
                val filteredNucleotides = fragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }

        val nucleotideCounts = SequenceCount.nucleotides(minBaseCount, filteredNucleotides)
        val nucleotideExonBoundaryStarts = aminoAcidBoundary.map { 3 * it }
        val homLoci = nucleotideCounts.homozygousIndices().toSet()
        val homStarts = nucleotideExonBoundaryStarts.filter { homLoci.contains(it) }
        val homEnds = nucleotideExonBoundaryStarts.filter { homLoci.contains(it + 1) && homLoci.contains(it + 2) }

        val result = mutableListOf<NucleotideFragment>()

        for (fragment in fragments) {
            var enriched = fragment
            for (homStart in homStarts) {
                if (missingStart(homStart, enriched)) {
                    enriched = enriched.addStart(homStart, nucleotideCounts)
                }
            }

            for (homEnd in homEnds) {
                if (missingEnd(homEnd, enriched)) {
                    enriched = enriched.addEnd(homEnd, nucleotideCounts)
                }
            }
            result.add(enriched)
        }

        return result

         */

        /*
        void $receiver$iv$iv;
        void $receiver$iv$iv2;
        void $receiver$iv$iv3;
        Object object;
        Object it;
        Collection collection;
        Iterable $receiver$iv$iv4;
        Iterable $receiver$iv;
        Intrinsics.checkParameterIsNotNull(fragments, (String) "fragments");
        Iterable iterable = $receiver$iv = (Iterable) fragments;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv4)
        {
            NucleotideFragment nucleotideFragment = (NucleotideFragment) item$iv$iv;
            collection = destination$iv$iv;
            boolean bl = false;
            object = ((NucleotideFragment) it).qualityFilter(mMinBaseQuality);
            collection.add(object);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv4 = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv4)
        {
            it = (NucleotideFragment) element$iv$iv;
            boolean bl = false;
            if(!((NucleotideFragment) it).isNotEmpty())
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List filteredNucleotides = (List) destination$iv$iv;
        SequenceCount nucleotideCounts = SequenceCount.Companion.nucleotides(mMinBaseCount, filteredNucleotides);
        Iterable $receiver$iv2 = mAminoAcidBoundary;
        Iterable iterable2 = $receiver$iv2;
        Collection destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv2, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv3)
        {
            void it2;
            int $i$f$filterTo = ((Number) item$iv$iv).intValue();
            collection = destination$iv$iv2;
            boolean bl = false;
            object = 3 * it2;
            collection.add(object);
        }
        List nucleotideExonBoundaryStarts = (List) destination$iv$iv2;
        Set homLoci = CollectionsKt.toSet((Iterable) nucleotideCounts.homozygousIndices());
        Iterable $receiver$iv3 = nucleotideExonBoundaryStarts;
        it = $receiver$iv3;
        Iterable destination$iv$iv3 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv2)
        {
            int it3 = ((Number) element$iv$iv).intValue();
            boolean bl = false;
            if(!homLoci.contains(it3))
            {
                continue;
            }
            destination$iv$iv3.add(element$iv$iv);
        }
        List homStarts = (List) destination$iv$iv3;
        Iterable $receiver$iv4 = nucleotideExonBoundaryStarts;
        destination$iv$iv3 = $receiver$iv4;
        Collection destination$iv$iv4 = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            int it4 = ((Number) element$iv$iv).intValue();
            boolean bl = false;
            if(!(homLoci.contains(it4 + 1) && homLoci.contains(it4 + 2)))
            {
                continue;
            }
            destination$iv$iv4.add(element$iv$iv);
        }
        List homEnds = (List) destination$iv$iv4;
        List result = new ArrayList();
        Iterator<? extends NucleotideFragment> iterator = fragments.iterator();
        while(iterator.hasNext())
        {
            NucleotideFragment fragment;
            NucleotideFragment enriched = fragment = iterator.next();
            Iterator iterator2 = homStarts.iterator();
            while(iterator2.hasNext())
            {
                int homStart = ((Number) iterator2.next()).intValue();
                if(!missingStart(homStart, enriched))
                {
                    continue;
                }
                enriched = addStart(enriched, homStart, nucleotideCounts);
            }
            iterator2 = homEnds.iterator();
            while(iterator2.hasNext())
            {
                int homEnd = ((Number) iterator2.next()).intValue();
                if(!missingEnd(homEnd, enriched))
                {
                    continue;
                }
                enriched = addEnd(enriched, homEnd, nucleotideCounts);
            }
            result.add(enriched);
        }
        return result;
         */
    }

    private boolean missingStart(int index, NucleotideFragment fragment)
    {
        return !fragment.containsNucleotide(index) && fragment.containsAllNucleotides(Sets.newHashSet(index + 1, index + 2));
    }

    private boolean missingEnd(int index, NucleotideFragment fragment)
    {
        return fragment.containsNucleotide(index) && !fragment.containsNucleotide(index + 1) && !fragment.containsNucleotide(index + 2);
    }

    private NucleotideFragment addStart(final NucleotideFragment fragment, int index, SequenceCount nucleotideCounts)
    {
        return fragment.enrich(index, nucleotideCounts.sequenceAt(index).get(0), mMinBaseQuality);
    }

    private NucleotideFragment addEnd(final NucleotideFragment fragment, int index, SequenceCount nucleotideCounts)
    {
        return fragment.enrich(
                index + 1, nucleotideCounts.sequenceAt(index + 1).get(0), mMinBaseQuality)
                .enrich(index + 2, nucleotideCounts.sequenceAt(index + 2).get(0), mMinBaseQuality);
    }

}
