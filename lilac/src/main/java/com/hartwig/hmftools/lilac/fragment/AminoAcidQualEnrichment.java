package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.MIN_DEPTH_FILTER;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline.RAW_REF_AMINO_ACID_COUNTS;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.copyNucleotideFragment;

import java.util.List;
import java.util.NavigableSet;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.evidence.AminoAcid;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public final class AminoAcidQualEnrichment
{
    public static List<Fragment> qualityFilterAminoAcidFragments(
            final HlaContext context, final List<Fragment> enrichedFragments, double minEvidenceFactor)
    {
        // only permit high quality amino acids, ie, amino acids that have at least [minEvidence]
        List<Fragment> qualityFilteredAminoAcidFragments = enrichedFragments.stream()
                .map(x -> copyNucleotideFragment(x))
                .collect(Collectors.toList());

        qualityFilteredAminoAcidFragments.forEach(x -> x.buildAminoAcids());

        SequenceCount highQualityAminoAcidCounts = SequenceCount.buildFromAminoAcids(minEvidenceFactor, qualityFilteredAminoAcidFragments);

        List<Fragment> unfilteredAminoAcidFragments = enrichedFragments.stream()
                .filter(x -> x.hasNucleotides())
                .map(x -> copyNucleotideFragment(x))
                .collect(Collectors.toList());

        unfilteredAminoAcidFragments.forEach(x -> x.buildAminoAcids());
        unfilteredAminoAcidFragments.forEach(x -> applyQualFilter(x, highQualityAminoAcidCounts, RAW_REF_AMINO_ACID_COUNTS.get(context.geneName())));

        return unfilteredAminoAcidFragments;
    }

    @VisibleForTesting
    public static void applyQualFilter(final Fragment fragment, final SequenceCount count, final SequenceCount rawAminoAcidCounts)
    {
        NavigableSet<Integer> filteredIntersect = Sets.newTreeSet();
        fragment.aminoAcidsByLoci().values().stream()
                .mapToInt(AminoAcid::locus)
                .forEach(locus ->
                {
                    if(rawAminoAcidCounts.get(locus).size() < MIN_DEPTH_FILTER)
                    {
                        filteredIntersect.add(locus);
                        return;
                    }

                    Set<String> allowed = Sets.newHashSet(count.getMinEvidenceSequences(locus));
                    String actual = fragment.aminoAcid(locus);
                    if(allowed.contains(actual))
                        filteredIntersect.add(locus);
                });

        fragment.filterAminoAcidsOnLoci(filteredIntersect);
    }
}
