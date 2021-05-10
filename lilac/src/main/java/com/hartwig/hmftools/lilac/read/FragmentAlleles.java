package com.hartwig.hmftools.lilac.read;

import static com.hartwig.hmftools.lilac.seq.HlaSequenceMatch.FULL;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceMatch.PARTIAL;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceMatch.WILD;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.HlaSequenceMatch;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class FragmentAlleles
{
    private final AminoAcidFragment mFragment;
    private final List<HlaAllele> mFull;
    private final List<HlaAllele> mPartial;
    private final List<HlaAllele> mWild;

    public FragmentAlleles(
            final AminoAcidFragment fragment, final List<HlaAllele> full, final List<HlaAllele> partial, final List<HlaAllele> wild)
    {
        mFragment = fragment;
        mFull = full;
        mPartial = partial;
        mWild = wild;
    }

    public final boolean contains(final HlaAllele allele)
    {
        return mFull.contains(allele) || mPartial.contains(allele) || mWild.contains(allele);
    }

    public final AminoAcidFragment getFragment() { return mFragment; }

    public final List<HlaAllele> getFull() { return mFull; }
    public final List<HlaAllele> getPartial() { return mPartial; }
    public final List<HlaAllele> getWild() { return mWild; }

    private static FragmentAlleles filter(final FragmentAlleles fragAlleles, final List<HlaAllele> alleles)
    {
        return new FragmentAlleles(
                fragAlleles.getFragment(),
                fragAlleles.getFull().stream().filter(x -> alleles.contains(x)).collect(Collectors.toList()),
                fragAlleles.getPartial().stream().filter(x -> alleles.contains(x)).collect(Collectors.toList()),
                fragAlleles.getWild().stream().filter(x -> alleles.contains(x)).collect(Collectors.toList()));
    }

    public static List<FragmentAlleles> filter(final List<FragmentAlleles> fragAlleleList, final List<HlaAllele> alleles)
    {
        return fragAlleleList.stream()
                .map(x -> x.filter(x, alleles))
                .filter(x -> !x.getFull().isEmpty() || !x.getPartial().isEmpty())
                .collect(Collectors.toList());
    }

    public static List<FragmentAlleles> create(
            final List<AminoAcidFragment> aminoAcidFragments, final List<Integer> hetLoci,
            final List<HlaSequenceLoci> sequences, final List<Integer> nucleotideLoci,
            final List<HlaSequenceLoci> nucleotideSequences)
    {
        return aminoAcidFragments.stream()
                .map(x -> create(x, hetLoci, sequences, nucleotideLoci, nucleotideSequences))
                .filter(x -> !x.getFull().isEmpty() || !x.getPartial().isEmpty())
                .collect(Collectors.toList());
    }

    private static FragmentAlleles create(
            final AminoAcidFragment aminoAcidFragment, final List<Integer> aminoAcidLoci, final List<HlaSequenceLoci> aminoAcidSequences,
            final List<Integer> nucleotideLoci, final List<HlaSequenceLoci> nucleotideSequences)
    {
        List<Integer> fragmentNucleotideLoci = aminoAcidFragment.getAminoAcidLoci().stream()
                .filter(x -> nucleotideLoci.contains(x)).collect(Collectors.toList());

        Collections.sort(fragmentNucleotideLoci);

        String fragmentNucleotides = aminoAcidFragment.nucleotides(fragmentNucleotideLoci);

        Map<HlaAllele,HlaSequenceMatch> alleleMatches = Maps.newHashMap();
        for(HlaSequenceLoci sequence : nucleotideSequences)
        {
            HlaAllele allele = sequence.getAllele();
            HlaSequenceMatch matchType = sequence.match(fragmentNucleotides, fragmentNucleotideLoci);
            if(matchType == HlaSequenceMatch.NONE)
                continue;

            // CHECK
            String gene = "HLA-" + allele.Gene;
            if(!aminoAcidFragment.getGenes().contains(gene))
                continue;

            alleleMatches.put(allele.asFourDigit(), matchType);
        }

        // CHECK prioritises the same way

        /*
        val matchingNucleotideSequences = nucleotideSequences
                .map { Pair(it.allele, it.match(fragmentNucleotides, *fragmentNucleotideLoci)) }
                    .filter { it.second != HlaSequenceMatch.NONE }
                    .filter { aminoAcidFragment.genes.contains("HLA-${it.first.gene}") }
                    .map { Pair(it.first.asFourDigit(), it.second) }
                    .distinct()
         */

        List<HlaAllele> fullNucleotideMatch = alleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> partialNucleotideMatch = alleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == PARTIAL || x.getValue() == WILD)
                .filter(x -> !fullNucleotideMatch.contains(x))
                .map(x -> x.getKey()).collect(Collectors.toList());

        // CHECK is sorting really required
        // TODO - should be intersection not union
        // val fragmentAminoAcidLoci = (aminoAcidFragment.aminoAcidLoci() intersect aminoAcidLoci).sorted().toIntArray()
        Set<Integer> combinedAminoAcidLoci = aminoAcidFragment.getAminoAcidLoci().stream().collect(Collectors.toSet());
        combinedAminoAcidLoci.addAll(aminoAcidLoci);
        List<Integer> fragmentAminoAcidLoci = combinedAminoAcidLoci.stream().collect(Collectors.toList());
        Collections.sort(fragmentAminoAcidLoci);

        String fragmentAminoAcids = aminoAcidFragment.aminoAcids(fragmentAminoAcidLoci);

        /*
        val matchingAminoAcidSequences = aminoAcidSequences.stream().map(x -> )
                .map { Pair(it.allele, it.match(fragmentAminoAcids, *fragmentAminoAcidLoci)) }
                    .filter { it.second != HlaSequenceMatch.NONE }
                    .filter { aminoAcidFragment.genes.contains("HLA-${it.first.gene}") }
         */

        List<HlaAllele> fullAminoAcidMatch = Lists.newArrayList();
        List<HlaAllele> partialAminoAcidMatch = Lists.newArrayList();
        List<HlaAllele> wildAminoAcidMatch = Lists.newArrayList();
        for(HlaSequenceLoci sequence : aminoAcidSequences)
        {
            HlaAllele allele = sequence.getAllele();
            HlaSequenceMatch matchType = sequence.match(fragmentAminoAcids, fragmentAminoAcidLoci);
            if(matchType == HlaSequenceMatch.NONE)
                continue;

            // CHECK
            String gene = "HLA-" + allele.Gene;
            if(!aminoAcidFragment.getGenes().contains(gene))
                continue;

            if(matchType == FULL)
                fullAminoAcidMatch.add(allele);
            else if(matchType == PARTIAL)
                partialAminoAcidMatch.add(allele);
            if(matchType == WILD)
                wildAminoAcidMatch.add(allele);
        }

        if (fullNucleotideMatch.isEmpty() && partialNucleotideMatch.isEmpty())
            return new FragmentAlleles(aminoAcidFragment, fullAminoAcidMatch, partialAminoAcidMatch, wildAminoAcidMatch);

        List<HlaAllele> consistentFull = fullAminoAcidMatch.stream()
                .filter(x -> fullNucleotideMatch.contains(x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> downgradedToPartial = fullAminoAcidMatch.stream()
                .filter(x -> partialNucleotideMatch.contains(x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> otherPartial = partialAminoAcidMatch.stream()
                .filter(x -> partialNucleotideMatch.contains(x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> distinctPartial = downgradedToPartial;
        otherPartial.stream().filter(x -> !downgradedToPartial.contains(x)).forEach(x -> distinctPartial.add(x));

        return new FragmentAlleles(aminoAcidFragment, consistentFull, distinctPartial, wildAminoAcidMatch);
    }
}
