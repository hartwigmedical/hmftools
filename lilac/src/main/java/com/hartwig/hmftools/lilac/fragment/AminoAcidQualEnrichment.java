package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.copyNucleotideFragment;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.lilac.seq.SequenceCount;

public final class AminoAcidQualEnrichment
{
    public static List<Fragment> qualityFilterAminoAcidFragments(final List<Fragment> fragments, double minEvidence)
    {
        // only permit high quality amino acids, ie, amino acids that have at least [minEvidence]
        List<Fragment> qualityFilteredAminoAcidFragments = fragments.stream()
                .map(x -> copyNucleotideFragment(x)).collect(Collectors.toList());

        qualityFilteredAminoAcidFragments.forEach(x -> x.buildAminoAcids());

        SequenceCount highQualityAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, qualityFilteredAminoAcidFragments);

        List<Fragment> unfilteredAminoAcidFragments = fragments.stream()
                .filter(x -> x.hasNucleotides())
                .map(x -> copyNucleotideFragment(x)).collect(Collectors.toList());

        unfilteredAminoAcidFragments.forEach(x -> x.buildAminoAcids());

        unfilteredAminoAcidFragments.forEach(x -> applyQualFilter(x, highQualityAminoAcidCounts));

        return unfilteredAminoAcidFragments;
    }

    private static void applyQualFilter(final Fragment fragment, final SequenceCount count)
    {
        List<Integer> initialIntersect = fragment.aminoAcidLoci();

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
