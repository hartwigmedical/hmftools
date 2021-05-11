package com.hartwig.hmftools.lilac.fragment;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.lilac.SequenceCount;

public class AminoAcidQualEnrichment
{
    private final int mMinEvidence;

    public AminoAcidQualEnrichment(int minEvidence)
    {
        mMinEvidence = minEvidence;
    }

    // Only permit high quality amino acids, ie, amino acids that have at least [minEvidence]
    public final List<AminoAcidFragment> enrich(final List<NucleotideFragment> fragments)
    {
        final List<AminoAcidFragment> qualityFilteredAminoAcidFragments
                = fragments.stream().map(x -> x.toAminoAcidFragment()).collect(Collectors.toList());

        SequenceCount highQualityAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, qualityFilteredAminoAcidFragments);

        final List<AminoAcidFragment> unfilteredAminoAcidFragments = fragments.stream()
                .filter(x -> x.isNotEmpty())
                .map(x -> x.toAminoAcidFragment())
                .collect(Collectors.toList());

        return unfilteredAminoAcidFragments.stream().map(x -> enrich(x, highQualityAminoAcidCounts)).collect(Collectors.toList());
    }

    private final AminoAcidFragment enrich(final AminoAcidFragment fragment, final SequenceCount count)
    {
        List<Integer> initialIntersect = fragment.getAminoAcidLoci();

        List<Integer> filteredIntersect = initialIntersect.stream()
                .filter(x -> filter(fragment, count, x))
                .collect(Collectors.toList());

        return fragment.intersectAminoAcidLoci(filteredIntersect);
    }

    private static boolean filter(final AminoAcidFragment fragment, final SequenceCount count, int loci)
    {
        List<String> allowed = count.getMinCountSequences(loci);
        String actual = fragment.aminoAcid(loci);

        return allowed.contains(actual);
    }

}
