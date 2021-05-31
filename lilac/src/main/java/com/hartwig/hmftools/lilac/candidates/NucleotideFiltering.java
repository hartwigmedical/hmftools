package com.hartwig.hmftools.lilac.candidates;

import static com.hartwig.hmftools.lilac.LilacConstants.GENE_IDS;
import static com.hartwig.hmftools.lilac.LilacConstants.getAminoAcidExonBoundaries;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class NucleotideFiltering
{
    private final int mMinNucleotideCount;
    private final List<Integer> mAminoAcidBoundaries;

    public NucleotideFiltering(int minNucleotideCount, final List<Integer> aminoAcidBoundaries)
    {
        mMinNucleotideCount = minNucleotideCount;
        mAminoAcidBoundaries = aminoAcidBoundaries;
    }

    public List<HlaSequenceLoci> filterCandidatesOnAminoAcidBoundaries(
            final List<HlaSequenceLoci> candidates, final List<Fragment> fragments)
    {
        List<HlaSequenceLoci> results = Lists.newArrayList();
        results.addAll(candidates);

        for(int boundary : mAminoAcidBoundaries)
        {
            int nucleotideStart = boundary * 3;
            final List<String> startSequences = nucleotideSequence(fragments, Lists.newArrayList(nucleotideStart));
            final List<String> endSequences = nucleotideSequence(fragments, Lists.newArrayList(nucleotideStart + 1, nucleotideStart + 2));

            results = results.stream()
                    .filter(x -> consistentWithAny(x, nucleotideStart, startSequences, endSequences))
                    .collect(Collectors.toList());
        }

        return results;
    }

    private static boolean consistentWithAny(
            final HlaSequenceLoci seqLoci, int startLoci, final List<String> startSequences, final List<String> endSequences)
    {
        return seqLoci.consistentWithAny(startSequences, Lists.newArrayList(startLoci))
            && seqLoci.consistentWithAny(endSequences, Lists.newArrayList(startLoci + 1, startLoci + 2));
    }

    private final List<String> nucleotideSequence(final List<Fragment> fragments, final List<Integer> nucleotideIndices)
    {
        Map<String,Integer> sequenceCounts = Maps.newHashMap();

        for(Fragment fragment : fragments)
        {
            if(!fragment.containsAllNucleotides(nucleotideIndices))
                continue;

            String nucleotides = fragment.nucleotides(nucleotideIndices);

            Integer count = sequenceCounts.get(nucleotides);
            sequenceCounts.put(nucleotides, count != null ? count + 1 : 1);
        }

        return sequenceCounts.entrySet().stream()
                .filter(x -> x.getValue() >= mMinNucleotideCount)
                .map(x -> x.getKey())
                .collect(Collectors.toList());

    }

    public static Map<String,List<Integer>> calcNucleotideHeterogygousLoci(final List<Integer> refNucleotideHetLoci)
    {
        // convert from amino acid exon boundaries to nucleotides for each gene
        Map<String,List<Integer>> hetLociMap = Maps.newHashMap();

        for(String gene : GENE_IDS)
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
                    refNucleotideHetLoci.stream().filter(x -> nucleotideExonBoundaries.contains(x)).collect(Collectors.toList()));

        }

        return hetLociMap;
    }
}
