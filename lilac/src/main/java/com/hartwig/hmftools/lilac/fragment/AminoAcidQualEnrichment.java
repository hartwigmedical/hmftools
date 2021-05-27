package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.copyNucleotideFragment;

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
    public List<Fragment> enrich(final List<Fragment> fragments)
    {
        List<Fragment> qualityFilteredAminoAcidFragments = fragments.stream()
                .map(x -> copyNucleotideFragment(x)).collect(Collectors.toList());

        qualityFilteredAminoAcidFragments.forEach(x -> x.buildAminoAcids());

        SequenceCount highQualityAminoAcidCounts = SequenceCount.aminoAcids(mMinEvidence, qualityFilteredAminoAcidFragments);

        List<Fragment> unfilteredAminoAcidFragments = fragments.stream()
                .filter(x -> x.hasNucleotides())
                .map(x -> copyNucleotideFragment(x)).collect(Collectors.toList());

        unfilteredAminoAcidFragments.forEach(x -> x.buildAminoAcids());

        unfilteredAminoAcidFragments.forEach(x -> enrich(x, highQualityAminoAcidCounts));

        return unfilteredAminoAcidFragments;
    }

    private void enrich(final Fragment fragment, final SequenceCount count)
    {
        List<Integer> initialIntersect = fragment.getAminoAcidLoci();

        List<Integer> filteredIntersect = initialIntersect.stream()
                .filter(x -> filter(fragment, count, x))
                .collect(Collectors.toList());

        fragment.filterOnLoci(filteredIntersect);
    }

    private static boolean filter(final Fragment fragment, final SequenceCount count, int loci)
    {
        List<String> allowed = count.getMinCountSequences(loci);
        String actual = fragment.aminoAcid(loci);

        return allowed.contains(actual);
    }

}
