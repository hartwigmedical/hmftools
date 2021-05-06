package com.hartwig.hmftools.lilac.fragment;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.SequenceCount;

public final class NucleotideQualEnrichment
{
    private final int mMinBaseQuality;
    private final int mMinEvidence;

    // Only permit high quality nucleotide values, ie, nucleotides that have at least 1 high quality (> [minBaseQuality])
    // and at least [minEvidence] instances in aggregate

    public final List<NucleotideFragment> enrich(final List<NucleotideFragment> fragments)
    {
        final List<NucleotideFragment> highQualFragments = fragments.stream()
                .map(x -> x.qualityFilter(mMinBaseQuality)).collect(Collectors.toList());

        /*
                val highQualityFragments = nucleotideFragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        val highQualityCounts = SequenceCount.nucleotides(1, highQualityFragments)
        val rawCounts = SequenceCount.nucleotides(minEvidence, nucleotideFragments)
        val result = nucleotideFragments.map { enrich(it, highQualityCounts, rawCounts) }
        return result

         */

        // TODO
        return highQualFragments;
    }

    private final NucleotideFragment enrich(NucleotideFragment fragment, SequenceCount highQualityCount, SequenceCount rawCount)
    {
        final List<Integer> nucleotideLoci = Lists.newArrayList();
        final List<Integer> nucleotideQuality = Lists.newArrayList();
        final List<String> nucleotides = Lists.newArrayList();

        for(Integer lociIndex : fragment.getNucleotideLoci())
        {
            int currentQuality = fragment.getNucleotideQuality().get(lociIndex);
            String fragmentNucleotide = fragment.getNucleotides().get(lociIndex);
            Collection<String> highQualitySequences = highQualityCount.sequenceAt(lociIndex);
            Collection<String> rawSequences = rawCount.sequenceAt(lociIndex);

            // TODO
            Collection<String> allowedSequences = Collections.EMPTY_LIST; // rawSequences intersect highQualitySequences
            // Set<String> allowedSequences = rawSequences.in.intersect((Iterable) rawSequences, (Iterable) highQualitySequences);

            if(allowedSequences.contains(fragmentNucleotide))
            {
                nucleotideLoci.add(lociIndex);
                nucleotideQuality.add(currentQuality);
                nucleotides.add(fragmentNucleotide);
            }
        }

        return new NucleotideFragment(fragment.getId(), fragment.getGenes(), nucleotideLoci, nucleotideQuality, nucleotides);
    }

    public NucleotideQualEnrichment(int minBaseQuality, int minEvidence)
    {
        this.mMinBaseQuality = minBaseQuality;
        this.mMinEvidence = minEvidence;
    }
}
