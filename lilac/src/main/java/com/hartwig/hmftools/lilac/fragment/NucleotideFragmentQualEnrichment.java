package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_NUCLEOTIDE_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_NUCLEOTIDE_HIGH_QUAL_EVIDENCE_FACTOR;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public final class NucleotideFragmentQualEnrichment
{
    public static List<Fragment> qualityFilterFragments(
            double minEvidence, double minHighQualEvidence, final List<Fragment> fragments, final List<Fragment> highQualFrags)
    {
        // fragments are all in nucleotide-space

        // filter fragments so that each nucleotide has at least 1 base at or above the min-qual threshold, and
        // X fragments (minEvidence) at that base with any qual
        SequenceCount highQualCounts = SequenceCount.nucleotides(minHighQualEvidence, highQualFrags);
        SequenceCount rawCounts = SequenceCount.nucleotides(minEvidence, fragments);

        return fragments.stream().map(x -> applyQualityFilter(x, highQualCounts, rawCounts)).collect(Collectors.toList());
    }

    private static Fragment applyQualityFilter(final Fragment fragment, final SequenceCount highQualityCount, final SequenceCount rawCount)
    {
        // checks whether all nucleotides have qual above the required level - if so return this fragment unch, otherwise build a
        // new fragment just with these filtered loci
        final List<Integer> filteredIndices = Lists.newArrayList();
        boolean allPresent = true;

        for(int i = 0; i < fragment.nucleotideLoci().size(); ++i)
        {
            int locusIndex = fragment.nucleotideLoci().get(i);
            String fragmentNucleotide = fragment.nucleotides().get(i);

            List<String> highQualitySequences = highQualityCount.getMinCountSequences(locusIndex, DEFAULT_MIN_NUCLEOTIDE_HIGH_QUAL_EVIDENCE_FACTOR);
            List<String> rawSequences = rawCount.getMinCountSequences(locusIndex, DEFAULT_MIN_NUCLEOTIDE_EVIDENCE_FACTOR);
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
            filteredLoci.add(fragment.nucleotideLoci().get(index));
            filteredQuality.add(fragment.nucleotideQuality().get(index));
            filteredNucleotides.add(fragment.nucleotides().get(index));
        }

        Fragment newFragment = new Fragment(
                fragment.reads().get(0), fragment.readGene(), fragment.genes(), filteredLoci, filteredQuality, filteredNucleotides);

        newFragment.addReads(fragment);

        return newFragment;
    }
}
