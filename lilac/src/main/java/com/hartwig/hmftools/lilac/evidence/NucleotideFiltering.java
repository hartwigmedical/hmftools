package com.hartwig.hmftools.lilac.evidence;

import static com.hartwig.hmftools.lilac.LilacConstants.MIN_DEPTH_FILTER;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.ReferenceData.getAminoAcidExonBoundaries;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_DRB3;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_DRB4;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
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
            final HlaGene gene, final Collection<HlaSequenceLoci> candidates, final List<Fragment> fragments)
    {
        List<HlaSequenceLoci> results = Lists.newArrayList();
        results.addAll(candidates);

        for(int boundary : mAminoAcidBoundaries)
        {
            int nucleotideStart = boundary * 3;
            if(localSpanCount(gene, fragments, Lists.newArrayList(nucleotideStart, nucleotideStart + 1, nucleotideStart + 2))
                    < MIN_DEPTH_FILTER)
                continue;

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

    private static List<String> nucleotideSequence(final Iterable<Fragment> fragments, final Collection<Integer> nucleotideIndices)
    {
        Set<String> nucleotideSequences = Sets.newHashSet();
        for(Fragment fragment : fragments)
        {
            if(!fragment.containsAllNucleotideLoci(nucleotideIndices))
                continue;

            String nucleotides = fragment.nucleotides(nucleotideIndices);
            nucleotideSequences.add(nucleotides);
        }

        return Lists.newArrayList(nucleotideSequences);
    }

    private static int localSpanCount(final HlaGene gene, final Iterable<Fragment> fragments, final Collection<Integer> nucleotideIndices)
    {
        int count = 0;
        for(Fragment fragment : fragments)
        {
            if(gene != HLA_DRB3 && gene != HLA_DRB4 && fragment.readGene() != gene)
                continue;

            if(!fragment.containsAllNucleotideLoci(nucleotideIndices))
                continue;

            count++;
        }

        return count;
    }

    public static Map<HlaGene, List<Integer>> calcNucleotideHeterogygousLoci(final Collection<Integer> refNucleotideHetLoci)
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
                    refNucleotideHetLoci.stream().filter(nucleotideExonBoundaries::contains).collect(Collectors.toList()));

        }

        return hetLociMap;
    }
}
