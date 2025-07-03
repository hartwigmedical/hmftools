package com.hartwig.hmftools.lilac.evidence;

import static java.lang.Math.ceil;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.ReferenceData.getAminoAcidExonBoundaries;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class NucleotideFiltering
{
    private final int mMinFilterDepth;
    private final double mMinFactor;
    private final Iterable<Integer> mAminoAcidBoundaries;

    public NucleotideFiltering(final int minFilterDepth, final double minFactor, final Iterable<Integer> aminoAcidBoundaries)
    {
        mMinFilterDepth = minFilterDepth;
        mMinFactor = minFactor;
        mAminoAcidBoundaries = aminoAcidBoundaries;
    }

    public List<HlaSequenceLoci> filterCandidatesOnAminoAcidBoundaries(final Collection<HlaSequenceLoci> candidates,
            final Iterable<Fragment> fragments)
    {
        List<HlaSequenceLoci> results = Lists.newArrayList();
        results.addAll(candidates);

        for(int boundary : mAminoAcidBoundaries)
        {
            int nucleotideStart = boundary * 3;
            final List<String> startSequences = nucleotideSequence(fragments, Lists.newArrayList(nucleotideStart));
            final List<String> endSequences = nucleotideSequence(fragments,
                    Lists.newArrayList(nucleotideStart + 1, nucleotideStart + 2));

            results = results.stream()
                    .filter(x -> consistentWithAny(x, nucleotideStart, startSequences, endSequences))
                    .collect(Collectors.toList());
        }

        return results;
    }

    private static boolean consistentWithAny(
            final HlaSequenceLoci seqLoci, final int startLoci, final List<String> startSequences,
            final List<String> endSequences)
    {
        return seqLoci.consistentWithAny(startSequences, Lists.newArrayList(startLoci))
                && seqLoci.consistentWithAny(endSequences, Lists.newArrayList(startLoci + 1, startLoci + 2));
    }

    private List<String> nucleotideSequence(final Iterable<Fragment> fragments, final Collection<Integer> nucleotideIndices)
    {
        Multiset<String> sequenceCounts = HashMultiset.create();
        int coverage = 0;

        for(Fragment fragment : fragments)
        {
            if(!fragment.containsAllNucleotideLoci(nucleotideIndices))
            {
                continue;
            }

            String nucleotides = fragment.nucleotides(nucleotideIndices);
            sequenceCounts.add(nucleotides);
            coverage++;
        }

        int minCount = (int) ceil(mMinFactor * coverage);
        return sequenceCounts.entrySet().stream()
                .filter(x -> x.getCount() >= minCount)
                .map(Multiset.Entry::getElement)
                .collect(Collectors.toList());
    }

    public static Map<String, List<Integer>> calcNucleotideHeterogygousLoci(final Collection<Integer> refNucleotideHetLoci)
    {
        // convert from amino acid exon boundaries to nucleotides for each gene
        Map<String, List<Integer>> hetLociMap = Maps.newHashMap();

        for(String gene : GENE_CACHE.GeneIds)
        {
            List<Integer> aminoAcidExonBoundaries = getAminoAcidExonBoundaries(gene);

            List<Integer> nucleotideExonBoundaries = Lists.newArrayList();

            for(Integer boundary : aminoAcidExonBoundaries)
            {
                nucleotideExonBoundaries.add(boundary * 3);
                nucleotideExonBoundaries.add(boundary * 3 + 1);
                nucleotideExonBoundaries.add(boundary * 3 + 2);
            }

            hetLociMap.put(gene,
                    refNucleotideHetLoci.stream().filter(nucleotideExonBoundaries::contains)
                            .collect(Collectors.toList()));

        }

        return hetLociMap;
    }
}
