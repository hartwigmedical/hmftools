package com.hartwig.hmftools.lilac.fragment;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.SequenceCount;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

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
        final List<NucleotideFragment> filteredNucleotides = fragments.stream()
                .map(x -> x.qualityFilter(mMinBaseQuality))
                .filter(x -> x.isNotEmpty())
                .collect(Collectors.toList());

        SequenceCount nucleotideCounts = SequenceCount.nucleotides(mMinBaseCount, filteredNucleotides);
        Set<Integer> nucleotideExonBoundaryStarts = mAminoAcidBoundary.stream().map(x -> x * 3).collect(Collectors.toSet());
        List<Integer> homLoci = nucleotideCounts.homozygousIndices();

        Set<Integer> homStarts = nucleotideExonBoundaryStarts.stream().filter(x -> homLoci.contains(x)).collect(Collectors.toSet());

        Set<Integer> homEnds = nucleotideExonBoundaryStarts.stream()
                .filter(x -> homLoci.contains(x + 1) && homLoci.contains(x + 2)).collect(Collectors.toSet());

        final List<NucleotideFragment> results = Lists.newArrayList();

        for (NucleotideFragment fragment : fragments)
        {
            for (Integer homStart : homStarts)
            {
                if (missingStart(homStart, fragment))
                {
                    addStart(fragment, homStart, nucleotideCounts);
                }
            }

            for (Integer homEnd : homEnds)
            {
                if (missingEnd(homEnd, fragment))
                    addEnd(fragment, homEnd, nucleotideCounts);
            }

            results.add(fragment);
        }

        return results;
    }

    private boolean missingStart(int index, NucleotideFragment fragment)
    {
        return !fragment.containsNucleotide(index) && fragment.containsAllNucleotides(Lists.newArrayList(index + 1, index + 2));
    }

    private boolean missingEnd(int index, NucleotideFragment fragment)
    {
        return fragment.containsNucleotide(index) && !fragment.containsNucleotide(index + 1) && !fragment.containsNucleotide(index + 2);
    }

    private void addStart(final NucleotideFragment fragment, int index, SequenceCount nucleotideCounts)
    {
        fragment.enrich(index, nucleotideCounts.getMinCountSequences(index).get(0), mMinBaseQuality);
    }

    private void addEnd(final NucleotideFragment fragment, int index, SequenceCount nucleotideCounts)
    {
        fragment.enrich(index + 1, nucleotideCounts.getMinCountSequences(index + 1).get(0), mMinBaseQuality);
        fragment.enrich(index + 2, nucleotideCounts.getMinCountSequences(index + 2).get(0), mMinBaseQuality);
    }
}
