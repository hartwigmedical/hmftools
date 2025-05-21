package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.copyNucleotideFragment;

import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public final class AminoAcidQualEnrichment
{
    public static List<Fragment> qualityFilterAminoAcidFragments(final LilacConfig config, final List<Fragment> fragments, double minEvidence)
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

        unfilteredAminoAcidFragments.forEach(x -> applyQualFilter(config, x, highQualityAminoAcidCounts));

        return unfilteredAminoAcidFragments;
    }

    private static void applyQualFilter(final LilacConfig config, final Fragment fragment, final SequenceCount count)
    {
        List<Integer> initialIntersect = fragment.aminoAcidLoci();

        List<Integer> filteredIntersect = Lists.newArrayList();
        for(Integer locusIndex : initialIntersect)
        {
            List<String> allowed = count.getMinCountOrVafSequences(locusIndex, config.MinAminoAcidEvidenceFactor);
            String actual = fragment.aminoAcid(locusIndex);

            if(allowed.contains(actual))
                filteredIntersect.add(locusIndex);
        }

        fragment.filterOnLoci(filteredIntersect);
    }
}
