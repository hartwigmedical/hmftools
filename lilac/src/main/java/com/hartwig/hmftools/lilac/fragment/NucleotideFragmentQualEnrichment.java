package com.hartwig.hmftools.lilac.fragment;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.SequenceCount;

public final class NucleotideFragmentQualEnrichment
{
    public static List<Fragment> enrich(final int minEvidence, final List<Fragment> fragments, final List<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        // filter fragments so that each nucleotide has at least 1 base at or above the min-qual threshold, and one
        // X fragments (minEvidence) at that base with any qual
        SequenceCount highQualCounts = SequenceCount.nucleotides(1, highQualFrags);
        SequenceCount rawCounts = SequenceCount.nucleotides(minEvidence, fragments);

        return fragments.stream().map(x -> enrich(x, highQualCounts, rawCounts)).collect(Collectors.toList());
    }

    private static Fragment enrich(final Fragment fragment, final SequenceCount highQualityCount, final SequenceCount rawCount)
    {
        final List<Integer> filteredIndices = Lists.newArrayList();
        boolean allPresent = true;

        for(int i = 0; i < fragment.getNucleotideLoci().size(); ++i)
        {
            int lociIndex = fragment.getNucleotideLoci().get(i);
            String fragmentNucleotide = fragment.getNucleotides().get(i);
            List<String> highQualitySequences = highQualityCount.getMinCountSequences(lociIndex);
            List<String> rawSequences = rawCount.getMinCountSequences(lociIndex);
            List<String> allowedSequences = highQualitySequences.stream().filter(x -> rawSequences.contains(x)).collect(Collectors.toList());

            if(allowedSequences.contains(fragmentNucleotide))
            {
                filteredIndices.add(i);
            }
            else
            {
                allPresent = false;
            }
        }

        if(allPresent)
            return fragment;

        int filteredCount = filteredIndices.size();
        final List<Integer> filteredLoci = Lists.newArrayListWithExpectedSize(filteredCount);
        final List<Integer> filteredQuality = Lists.newArrayListWithExpectedSize(filteredCount);
        final List<String> filteredNucleotides = Lists.newArrayListWithExpectedSize(filteredCount);

        for(Integer index : filteredIndices)
        {
            filteredLoci.add(fragment.getNucleotideLoci().get(index));
            filteredQuality.add(fragment.getNucleotideQuality().get(index));
            filteredNucleotides.add(fragment.getNucleotides().get(index));
        }

        return new Fragment(fragment.id(), fragment.readInfo(), fragment.getGenes(), filteredLoci, filteredQuality, filteredNucleotides);
    }
}
