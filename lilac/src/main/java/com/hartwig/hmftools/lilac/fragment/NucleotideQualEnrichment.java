package com.hartwig.hmftools.lilac.fragment;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.SequenceCount;

public final class NucleotideQualEnrichment
{
    private final int mMinBaseQuality;
    private final int mMinEvidence;

    public NucleotideQualEnrichment(int minBaseQuality, int minEvidence)
    {
        mMinBaseQuality = minBaseQuality;
        mMinEvidence = minEvidence;
    }

    public final List<NucleotideFragment> enrich(final List<NucleotideFragment> fragments)
    {
        // Only permit high quality nucleotide values, ie, nucleotides that have at least 1 high quality (> [minBaseQuality])
        // and at least [minEvidence] instances in aggregate
        final List<NucleotideFragment> highQualFragments = fragments.stream()
                .map(x -> x.qualityFilter(mMinBaseQuality))
                .filter(x -> x.isNotEmpty())
                .collect(Collectors.toList());

        SequenceCount highQualCounts = SequenceCount.nucleotides(1, highQualFragments);
        SequenceCount rawCounts = SequenceCount.nucleotides(mMinEvidence, fragments);

        return fragments.stream().map(x -> enrich(x, highQualCounts, rawCounts)).collect(Collectors.toList());
    }

    private NucleotideFragment enrich(
            final NucleotideFragment fragment, final SequenceCount highQualityCount, final SequenceCount rawCount)
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

        return new NucleotideFragment(fragment.getId(), fragment.getGenes(), filteredLoci, filteredQuality, filteredNucleotides);
    }
}
