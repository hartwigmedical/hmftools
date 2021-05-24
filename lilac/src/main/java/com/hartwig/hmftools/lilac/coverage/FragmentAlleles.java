package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.seq.HlaSequenceMatch.FULL;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceMatch.PARTIAL;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceMatch.WILD;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.HlaSequenceMatch;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
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

    public boolean contains(final HlaAllele allele)
    {
        return mFull.contains(allele) || mPartial.contains(allele) || mWild.contains(allele);
    }

    public final AminoAcidFragment getFragment() { return mFragment; }

    public final List<HlaAllele> getFull() { return mFull; }
    public final List<HlaAllele> getPartial() { return mPartial; }
    public final List<HlaAllele> getWild() { return mWild; }

    public static List<FragmentAlleles> filter(final List<FragmentAlleles> fragAlleleList, final List<HlaAllele> alleles)
    {
        // gather any fragment allele which contains at least one of the specified alleles in its full or partial list,
        // then collecting any matching alleles in each of the three groups
        List<FragmentAlleles> matchedFragAlleles = Lists.newArrayList();

        for(FragmentAlleles fragAlleles : fragAlleleList)
        {
            FragmentAlleles matchedFragAllele = null;

            for(HlaAllele allele : alleles)
            {
                boolean inFull = fragAlleles.getFull().contains(allele);
                boolean inPartial = fragAlleles.getPartial().contains(allele);

                if((inFull || inPartial) && matchedFragAllele == null)
                {
                    matchedFragAllele = new FragmentAlleles(
                            fragAlleles.getFragment(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList());

                    matchedFragAlleles.add(matchedFragAllele);
                }

                if(matchedFragAllele != null)
                {
                    if(inFull)
                        matchedFragAllele.getFull().add(allele);

                    if(inPartial)
                        matchedFragAllele.getPartial().add(allele);
                }
            }

            if(matchedFragAllele != null)
            {
                List<HlaAllele> wildAlleles = alleles.stream().filter(x -> fragAlleles.getWild().contains(x)).collect(Collectors.toList());
                matchedFragAllele.getWild().addAll(wildAlleles);
            }
        }

        return matchedFragAlleles;
    }

    public static List<FragmentAlleles> createFragmentAlleles(
            final List<AminoAcidFragment> refCoverageFragments, final List<Integer> refAminoAcidHetLoci,
            final List<HlaSequenceLoci> candidateAminoAcidSequences, final List<Set<String>> refAminoAcids,
            final Map<String,List<Integer>> refNucleotideHetLoci, final List<HlaSequenceLoci> candidateNucleotideSequences,
            final List<Set<String>> refNucleotides)
    {
        List<FragmentAlleles> results = Lists.newArrayList();

        for(AminoAcidFragment fragment : refCoverageFragments)
        {
            FragmentAlleles fragmentAlleles = create(
                    fragment, refAminoAcidHetLoci, candidateAminoAcidSequences, refAminoAcids,
                    refNucleotideHetLoci, candidateNucleotideSequences, refNucleotides);

            if(!fragmentAlleles.getFull().isEmpty() || !fragmentAlleles.getPartial().isEmpty())
                results.add(fragmentAlleles);
        }

        return results;
    }

    private static FragmentAlleles create(
            final AminoAcidFragment fragment, final List<Integer> aminoAcidLoci, final List<HlaSequenceLoci> aminoAcidSequences,
            final List<Set<String>> refAminoAcids,
            final Map<String,List<Integer>> refNucleotideLoci, final List<HlaSequenceLoci> nucleotideSequences,
            final List<Set<String>> refNucleotides)
    {
        Map<String,List<Integer>> fragNucleotideLociMap = Maps.newHashMap();

        refNucleotideLoci.entrySet().forEach(x -> fragNucleotideLociMap.put(
                x.getKey(), fragment.getNucleotideLoci().stream().filter(y -> x.getValue().contains(y)).collect(Collectors.toList())));

        final List<HlaAllele> fullNucleotideMatch = Lists.newArrayList();
        final List<HlaAllele> partialNucleotideMatch = Lists.newArrayList();

        Map<HlaAllele,HlaSequenceMatch> nucleotideAlleleMatches = findNucleotideMatches(
                fragment, refNucleotideLoci, nucleotideSequences, refNucleotides);

        fullNucleotideMatch.addAll(nucleotideAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList()));

        partialNucleotideMatch.addAll(nucleotideAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == PARTIAL || x.getValue() == WILD)
                .filter(x -> !fullNucleotideMatch.contains(x))
                .map(x -> x.getKey()).collect(Collectors.toList()));

        Map<HlaAllele,HlaSequenceMatch> aminoAcidAlleleMatches = findAminoAcidMatches(
                fragment, aminoAcidLoci, aminoAcidSequences, refAminoAcids);

        List<HlaAllele> fullAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> partialAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == PARTIAL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> wildAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == WILD).map(x -> x.getKey()).collect(Collectors.toList());

        if(fullNucleotideMatch.isEmpty() && partialNucleotideMatch.isEmpty())
            return new FragmentAlleles(fragment, fullAminoAcidMatch, partialAminoAcidMatch, wildAminoAcidMatch);

        List<HlaAllele> consistentFull = fullAminoAcidMatch.stream()
                .filter(x -> fullNucleotideMatch.contains(x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> downgradedToPartial = fullAminoAcidMatch.stream()
                .filter(x -> partialNucleotideMatch.contains(x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> otherPartial = partialAminoAcidMatch.stream()
                .filter(x -> partialNucleotideMatch.contains(x.asFourDigit())).collect(Collectors.toList());

        List<HlaAllele> distinctPartial = downgradedToPartial;
        otherPartial.stream().filter(x -> !downgradedToPartial.contains(x)).forEach(x -> distinctPartial.add(x));

        return new FragmentAlleles(fragment, consistentFull, distinctPartial, wildAminoAcidMatch);
    }

    private static Map<HlaAllele,HlaSequenceMatch> findNucleotideMatches(
            final AminoAcidFragment fragment, final Map<String,List<Integer>> refNucleotideLociMap,
            final List<HlaSequenceLoci> nucleotideSequences, final List<Set<String>> refNucleotides)
    {
        Map<String,List<Integer>> fragNucleotideLociMap = Maps.newHashMap();

        //refNucleotideLociMap.entrySet().forEach(x -> fragNucleotideLociMap.put(
        //        x.getKey(), fragment.getNucleotideLoci().stream().filter(y -> x.getValue().contains(y)).collect(Collectors.toList())));

        Map<HlaAllele,HlaSequenceMatch> alleleMatches = Maps.newHashMap();

        // also attempt to retrive amino acids from low-qual nucleotides
        Map<Integer,String> missedNucleotides = Maps.newHashMap();

        for(Map.Entry<String,List<Integer>> entry : refNucleotideLociMap.entrySet())
        {
            String gene = entry.getKey();
            List<Integer> refNucleotideLoci = entry.getValue();

            List<Integer> fragmentMatchedLoci = fragment.getNucleotideLoci().stream()
                    .filter(y -> refNucleotideLoci.contains(y)).collect(Collectors.toList());

            // also check for support from low-qual reads at this same location
            List<Integer> missedNucleotideLoci = refNucleotideLoci.stream()
                    .filter(x -> !fragment.getNucleotideLoci().contains(x)).collect(Collectors.toList());

            for(Integer missedLocus : missedNucleotideLoci)
            {
                if(missedLocus >= refNucleotides.size())
                    continue;

                Set<String> candidateNucleotides = refNucleotides.get(missedLocus);

                String lowQualNucleotide = fragment.getLowQualNucleotide(missedLocus);

                if(!lowQualNucleotide.isEmpty() && candidateNucleotides.contains(lowQualNucleotide))
                {
                    fragmentMatchedLoci.add(missedLocus);
                    missedNucleotides.put(missedLocus, lowQualNucleotide);
                }
            }

            fragNucleotideLociMap.put(gene, fragmentMatchedLoci);
        }

        if(fragNucleotideLociMap.values().stream().allMatch(x -> x.isEmpty()))
            return alleleMatches;

        fragNucleotideLociMap.values().forEach(x -> Collections.sort(x));

        Map<String,String> fragNucleotideSequences = Maps.newHashMap();

        for(HlaSequenceLoci sequence : nucleotideSequences)
        {
            HlaAllele allele = sequence.Allele;
            if(!fragment.getGenes().contains(allele.geneName()))
                continue;

            List<Integer> fragNucleotideLoci = fragNucleotideLociMap.get(allele.Gene);
            if(fragNucleotideLoci.isEmpty())
                continue;

            String fragNucleotides = fragNucleotideSequences.get(allele.Gene);

            if(fragNucleotides == null)
            {
                if(missedNucleotides.isEmpty())
                {
                    fragNucleotides = fragment.nucleotides(fragNucleotideLoci);
                }
                else
                {
                    StringJoiner nucSj = new StringJoiner("");

                    for(Integer locus : fragNucleotideLoci)
                    {
                        if(missedNucleotides.containsKey(locus))
                            nucSj.add(missedNucleotides.get(locus));
                        else
                            nucSj.add(fragment.nucleotide(locus));
                    }

                    fragNucleotides = nucSj.toString();
                }

                fragNucleotideSequences.put(allele.Gene, fragNucleotides);
            }

            HlaSequenceMatch matchType = sequence.match(fragNucleotides, fragNucleotideLoci);
            if(matchType == HlaSequenceMatch.NONE)
                continue;

            // keep the best match
            Map.Entry<HlaAllele,HlaSequenceMatch> entryMatch = alleleMatches.entrySet().stream()
                    .filter(x -> x.getKey() == allele.asFourDigit()).findFirst().orElse(null);

            if(entryMatch == null || matchType.isBetter(entryMatch.getValue()))
            {
                if(entryMatch != null)
                    alleleMatches.remove(entryMatch.getKey());

                alleleMatches.put(allele.asFourDigit(), matchType);
            }
        }

        return alleleMatches;
    }

    private static Map<HlaAllele,HlaSequenceMatch> findAminoAcidMatches(
            final AminoAcidFragment fragment, final List<Integer> aminoAcidLoci, final List<HlaSequenceLoci> aminoAcidSequences,
            final List<Set<String>> refAminoAcids)
    {
        Map<HlaAllele,HlaSequenceMatch> alleleMatches = Maps.newHashMap();

        List<Integer> fragmentAminoAcidLoci = fragment.getAminoAcidLoci().stream()
                .filter(x -> aminoAcidLoci.contains(x)).collect(Collectors.toList());

        if(fragmentAminoAcidLoci.isEmpty())
            return alleleMatches;

        // also attempt to retrive amino acids from low-qual nucleotides
        Map<Integer,String> missedAminoAcids = Maps.newHashMap();

        List<Integer> missedAminoAcidLoci = aminoAcidLoci.stream()
                .filter(x -> !fragmentAminoAcidLoci.contains(x)).collect(Collectors.toList());

        for(Integer missedLocus : missedAminoAcidLoci)
        {
            if(missedLocus >= refAminoAcids.size())
                continue;

            Set<String> candidateAminoAcids = refAminoAcids.get(missedLocus);

            String lowQualAminoAcid = fragment.getLowQualAminoAcid(missedLocus);

            if(!lowQualAminoAcid.isEmpty() && candidateAminoAcids.contains(lowQualAminoAcid))
            {
                fragmentAminoAcidLoci.add(missedLocus);
                missedAminoAcids.put(missedLocus, lowQualAminoAcid);
            }
        }

        Collections.sort(fragmentAminoAcidLoci);

        String fragmentAminoAcids;

        if(missedAminoAcids.isEmpty())
        {
            fragmentAminoAcids = fragment.aminoAcids(fragmentAminoAcidLoci);
        }
        else
        {
            StringJoiner aaSj = new StringJoiner("");

            for(Integer locus : fragmentAminoAcidLoci)
            {
                if(missedAminoAcids.containsKey(locus))
                    aaSj.add(missedAminoAcids.get(locus));
                else
                    aaSj.add(fragment.aminoAcid(locus));
            }

            fragmentAminoAcids = aaSj.toString();
        }

        for(HlaSequenceLoci sequence : aminoAcidSequences)
        {
            HlaSequenceMatch matchType = sequence.match(fragmentAminoAcids, fragmentAminoAcidLoci);
            if(matchType == HlaSequenceMatch.NONE)
                continue;

            HlaAllele allele = sequence.Allele;

            if(!fragment.getGenes().contains(allele.geneName()))
                continue;

            Map.Entry<HlaAllele,HlaSequenceMatch> entryMatch = alleleMatches.entrySet().stream()
                    .filter(x -> x.getKey() == allele).findFirst().orElse(null);

            if(entryMatch == null || matchType.isBetter(entryMatch.getValue()))
            {
                if(entryMatch != null)
                    alleleMatches.remove(entryMatch.getKey());

                alleleMatches.put(allele, matchType);
            }
        }

        return alleleMatches;
    }

    public static void applyUniqueStopLossFragments(
            final List<FragmentAlleles> fragmentAlleles, int stopLossFragments, final List<HlaAllele> stopLossAlleles)
    {
        if(stopLossFragments == 0)
            return;

        for(HlaAllele stopLossAllele : stopLossAlleles)
        {
            List<FragmentAlleles> sampleFragments = fragmentAlleles.stream()
                    .filter(x -> x.contains(stopLossAllele))
                    .map(x -> new FragmentAlleles(x.getFragment(), Lists.newArrayList(stopLossAllele), Lists.newArrayList(), Lists.newArrayList()))
                    .collect(Collectors.toList());

            for(int i = 0; i < stopLossFragments; ++i)
            {
                if(i >= sampleFragments.size())
                    break;

                fragmentAlleles.add(sampleFragments.get(i));
            }
        }
    }
}
