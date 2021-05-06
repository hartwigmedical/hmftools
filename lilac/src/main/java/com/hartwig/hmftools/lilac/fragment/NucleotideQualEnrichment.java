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
                .map(x -> x.qualityFilter(mMinBaseQuality))
                .filter(x -> x.isNotEmpty())
                .collect(Collectors.toList());

        SequenceCount highQualCounts = SequenceCount.nucleotides(1, highQualFragments);
        SequenceCount rawCounts = SequenceCount.nucleotides(mMinEvidence, fragments);
        return fragments.stream().map(x -> enrich(x, highQualCounts, rawCounts)).collect(Collectors.toList());
    }

    private final NucleotideFragment enrich(
            final NucleotideFragment fragment, final SequenceCount highQualityCount, final SequenceCount rawCount)
    {
        final List<Integer> nucleotideLoci = Lists.newArrayList();
        final List<Integer> nucleotideQuality = Lists.newArrayList();
        final List<String> nucleotides = Lists.newArrayList();

        for(Integer lociIndex : fragment.getNucleotideLoci())
        {
            int currentQuality = fragment.getNucleotideQuality().get(lociIndex);
            String fragmentNucleotide = fragment.getNucleotides().get(lociIndex);
            List<String> highQualitySequences = highQualityCount.sequenceAt(lociIndex);
            List<String> rawSequences = rawCount.sequenceAt(lociIndex);
            List<String> allowedSequences = highQualitySequences.stream().filter(x -> rawSequences.contains(x)).collect(Collectors.toList());

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
