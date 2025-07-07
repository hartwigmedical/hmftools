package com.hartwig.hmftools.lilac.fragment;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.utils.AminoAcid;

public final class AminoAcidQualEnrichment
{
    private AminoAcidQualEnrichment() {}

    public static List<Fragment> qualityFilterAminoAcidFragments(final LilacConfig config, final Collection<Fragment> fragments)
    {
        // only permit high quality amino acids, ie, amino acids that have at least [minEvidence]
        List<Fragment> qualityFilteredAminoAcidFragments = fragments.stream()
                .map(FragmentUtils::copyNucleotideFragment).collect(Collectors.toList());

        qualityFilteredAminoAcidFragments.forEach(Fragment::buildAminoAcids);

        SequenceCount highQualityAminoAcidCounts = SequenceCount.aminoAcids(
		config.MinVafFilterDepth, config.MinEvidenceFactor, qualityFilteredAminoAcidFragments);

        List<Fragment> unfilteredAminoAcidFragments = fragments.stream()
                .filter(Fragment::hasNucleotides)
                .map(FragmentUtils::copyNucleotideFragment).collect(Collectors.toList());

        unfilteredAminoAcidFragments.forEach(Fragment::buildAminoAcids);

        unfilteredAminoAcidFragments.forEach(x -> applyQualFilter(config, x, highQualityAminoAcidCounts));

        return unfilteredAminoAcidFragments;
    }

    private static void applyQualFilter(final LilacConfig config, final Fragment fragment, final SequenceCount count)
    {
        Set<Integer> filteredIntersect = Sets.newHashSet();
        fragment.aminoAcidsByLoci().values().stream()
                .mapToInt(AminoAcid::locus)
                .forEach(locusIndex ->
                {
                    List<String> allowed = count.getMinEvidenceSequences(locusIndex, config.MinEvidenceFactor);
                    String actual = fragment.aminoAcid(locusIndex);

                    if(allowed.contains(actual))
                    {
                        filteredIntersect.add(locusIndex);
                    }
                });

        fragment.filterAminoAcidsOnLoci(filteredIntersect);
    }
}
