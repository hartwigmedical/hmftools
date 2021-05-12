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
                fragAlleles.getFull().stream().filter(x -> HlaAllele.contains(alleles, x)).collect(Collectors.toList()),
                fragAlleles.getPartial().stream().filter(x -> HlaAllele.contains(alleles, x)).collect(Collectors.toList()),
                fragAlleles.getWild().stream().filter(x -> HlaAllele.contains(alleles, x)).collect(Collectors.toList()));
    }

    public static List<FragmentAlleles> filter(final List<FragmentAlleles> fragAlleleList, final List<HlaAllele> alleles)
    {
        return fragAlleleList.stream()
                .map(x -> x.filter(x, alleles))
                .filter(x -> !x.getFull().isEmpty() || !x.getPartial().isEmpty())
                .collect(Collectors.toList());
    }

    public static List<FragmentAlleles> create(
            final List<AminoAcidFragment> referenceCoverageFragments, final List<Integer> referenceAminoAcidHeterozygousLoci,
            final List<HlaSequenceLoci> candidateAminoAcidSequences, final List<Integer> referenceNucleotideHeterozygousLoci,
            final List<HlaSequenceLoci> candidateNucleotideSequences)
    {
        List<FragmentAlleles> results = Lists.newArrayList();

        for(AminoAcidFragment fragment : referenceCoverageFragments)
        {
            FragmentAlleles fragmentAlleles = create(
                    fragment, referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences, referenceNucleotideHeterozygousLoci, candidateNucleotideSequences);

            if(!fragmentAlleles.getFull().isEmpty() || !fragmentAlleles.getPartial().isEmpty())
                results.add(fragmentAlleles);
        }

        return results;
    }

    private static FragmentAlleles create(
            final AminoAcidFragment aminoAcidFragment, final List<Integer> aminoAcidLoci, final List<HlaSequenceLoci> aminoAcidSequences,
            final List<Integer> nucleotideLoci, final List<HlaSequenceLoci> nucleotideSequences)
    {
        List<Integer> fragmentNucleotideLoci = aminoAcidFragment.getNucleotideLoci().stream()
                .filter(x -> nucleotideLoci.contains(x)).collect(Collectors.toList());

        Map<HlaAllele,HlaSequenceMatch> alleleMatches = Maps.newHashMap();

        if(!fragmentNucleotideLoci.isEmpty())
        {
            Collections.sort(fragmentNucleotideLoci);

            String fragmentNucleotides = aminoAcidFragment.nucleotides(fragmentNucleotideLoci);

            for(HlaSequenceLoci sequence : nucleotideSequences)
            {
                HlaSequenceMatch matchType = sequence.match(fragmentNucleotides, fragmentNucleotideLoci);
                if(matchType == HlaSequenceMatch.NONE)
                    continue;

                HlaAllele allele = sequence.getAllele();
                String gene = "HLA-" + allele.Gene;
                if(!aminoAcidFragment.getGenes().contains(gene))
                    continue;

                // keep the best match
                Map.Entry<HlaAllele,HlaSequenceMatch> entryMatch = alleleMatches.entrySet().stream()
                        .filter(x -> x.getKey().matches(allele.asFourDigit())).findFirst().orElse(null);

                if(entryMatch == null || matchType.isBetter(entryMatch.getValue()))
                {
                    if(entryMatch != null)
                        alleleMatches.remove(entryMatch.getKey());

                    alleleMatches.put(allele.asFourDigit(), matchType);
                }
            }
        }

        List<HlaAllele> fullNucleotideMatch = alleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> partialNucleotideMatch = alleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == PARTIAL || x.getValue() == WILD)
                .filter(x -> !fullNucleotideMatch.contains(x))
                .map(x -> x.getKey()).collect(Collectors.toList());

        List<Integer> fragmentAminoAcidLoci = aminoAcidFragment.getAminoAcidLoci().stream()
                .filter(x -> aminoAcidLoci.contains(x)).collect(Collectors.toList());

        alleleMatches.clear();

        if(!fragmentAminoAcidLoci.isEmpty())
        {
            Collections.sort(fragmentAminoAcidLoci);

            String fragmentAminoAcids = aminoAcidFragment.aminoAcids(fragmentAminoAcidLoci);

            /*
            val matchingAminoAcidSequences = aminoAcidSequences.stream().map(x -> )
                    .map { Pair(it.allele, it.match(fragmentAminoAcids, *fragmentAminoAcidLoci)) }
                        .filter { it.second != HlaSequenceMatch.NONE }
                        .filter { aminoAcidFragment.genes.contains("HLA-${it.first.gene}") }
             */

            for(HlaSequenceLoci sequence : aminoAcidSequences)
            {
                HlaSequenceMatch matchType = sequence.match(fragmentAminoAcids, fragmentAminoAcidLoci);
                if(matchType == HlaSequenceMatch.NONE)
                    continue;

                HlaAllele allele = sequence.getAllele();

                String gene = "HLA-" + allele.Gene;
                if(!aminoAcidFragment.getGenes().contains(gene))
                    continue;

                Map.Entry<HlaAllele,HlaSequenceMatch> entryMatch = alleleMatches.entrySet().stream()
                        .filter(x -> x.getKey().matches(allele)).findFirst().orElse(null);

                if(entryMatch == null || matchType.isBetter(entryMatch.getValue()))
                {
                    if(entryMatch != null)
                        alleleMatches.remove(entryMatch.getKey());

                    alleleMatches.put(allele, matchType);
                }
            }
        }

        List<HlaAllele> fullAminoAcidMatch = alleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> partialAminoAcidMatch = alleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == PARTIAL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> wildAminoAcidMatch = alleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == WILD).map(x -> x.getKey()).collect(Collectors.toList());

        if (fullNucleotideMatch.isEmpty() && partialNucleotideMatch.isEmpty())
            return new FragmentAlleles(aminoAcidFragment, fullAminoAcidMatch, partialAminoAcidMatch, wildAminoAcidMatch);

        List<HlaAllele> consistentFull = fullAminoAcidMatch.stream()
                .filter(x -> HlaAllele.contains(fullNucleotideMatch, x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> downgradedToPartial = fullAminoAcidMatch.stream()
                .filter(x -> HlaAllele.contains(partialNucleotideMatch, x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> otherPartial = partialAminoAcidMatch.stream()
                .filter(x -> HlaAllele.contains(partialNucleotideMatch, x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> distinctPartial = downgradedToPartial;
        otherPartial.stream().filter(x -> !HlaAllele.contains(downgradedToPartial, x)).forEach(x -> distinctPartial.add(x));

        return new FragmentAlleles(aminoAcidFragment, consistentFull, distinctPartial, wildAminoAcidMatch);
    }
}
