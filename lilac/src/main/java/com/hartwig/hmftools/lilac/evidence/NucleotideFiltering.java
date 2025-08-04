package com.hartwig.hmftools.lilac.evidence;

import static java.lang.Math.ceil;
import static java.lang.Math.max;

import static com.hartwig.hmftools.lilac.LilacConstants.MIN_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_EVIDENCE_SUPPORT;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.ReferenceData.getAminoAcidExonBoundaries;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

public class NucleotideFiltering
{
    private final List<Integer> mAminoAcidBoundaries;

    public NucleotideFiltering(final List<Integer> aminoAcidBoundaries)
    {
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

    private List<String> nucleotideSequence(final List<Fragment> fragments, final List<Integer> nucleotideIndices)
    {
        Map<String, Integer> sequenceCounts = Maps.newHashMap();

        int totalCount = 0;
        for(Fragment fragment : fragments)
        {
            if(!fragment.containsAllNucleotideLoci(nucleotideIndices))
                continue;

            String nucleotides = fragment.nucleotides(nucleotideIndices);
            sequenceCounts.merge(nucleotides, 1, Integer::sum);
            totalCount++;
        }

        int minNucleotideCount = max(MIN_EVIDENCE_SUPPORT, (int)ceil(totalCount * MIN_EVIDENCE_FACTOR));

        return sequenceCounts.entrySet().stream()
                .filter(x -> x.getValue() >= minNucleotideCount)
                .map(x -> x.getKey())
                .collect(Collectors.toList());

    }

    public static Map<HlaGene, List<Integer>> calcNucleotideHeterogygousLoci(final List<Integer> refNucleotideHetLoci)
    {
        // convert from amino acid exon boundaries to nucleotides for each gene
        Map<HlaGene, List<Integer>> hetLociMap = Maps.newHashMap();

        for(HlaGene gene : GENE_CACHE.GeneNames)
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
