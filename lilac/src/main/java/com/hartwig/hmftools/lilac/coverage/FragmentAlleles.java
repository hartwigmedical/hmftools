package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_Y_FRAGMENT_THRESHOLD;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.filterExonBoundaryWildcards;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.filterWildcards;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.FULL;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.WILD;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class FragmentAlleles
{
    private final Fragment mFragment;
    private final List<HlaAllele> mFull;
    private final List<HlaAllele> mWild;

    public FragmentAlleles(
            final Fragment fragment, final List<HlaAllele> full, final List<HlaAllele> wild)
    {
        mFragment = fragment;
        mFull = full;
        mWild = wild;
    }

    public boolean contains(final HlaAllele allele)
    {
        return mFull.contains(allele) || mWild.contains(allele);
    }

    public final Fragment getFragment() { return mFragment; }

    public final List<HlaAllele> getFull() { return mFull; }
    public final List<HlaAllele> getWild() { return mWild; }

    public static List<FragmentAlleles> filter(final List<FragmentAlleles> fragAlleleList, final List<HlaAllele> alleles)
    {
        // gather any fragment allele which contains at least one of the specified alleles in its full or wild list,
        // then collecting any matching alleles in each of the three groups
        List<FragmentAlleles> matchedFragAlleles = Lists.newArrayList();

        for(FragmentAlleles fragAllele : fragAlleleList)
        {
            if(alleles.stream().anyMatch(x -> fragAllele.contains(x)))
            {
                matchedFragAlleles.add(new FragmentAlleles(
                        fragAllele.getFragment(),
                        alleles.stream().filter(x -> fragAllele.getFull().contains(x)).collect(Collectors.toList()),
                        alleles.stream().filter(x -> fragAllele.getWild().contains(x)).collect(Collectors.toList())));
            }
        }

        return matchedFragAlleles;
    }

    public static List<FragmentAlleles> createFragmentAlleles(
            final List<Fragment> refCoverageFragments, final List<Integer> refAminoAcidHetLoci,
            final List<HlaSequenceLoci> candidateAminoAcidSequences, final List<Set<String>> refAminoAcids,
            final Map<String,List<Integer>> refNucleotideHetLoci, final List<HlaSequenceLoci> candidateNucleotideSequences,
            final List<Set<String>> refNucleotides)
    {
        LL_LOGGER.info("building frag-alleles from aminoAcids(frags={} hetLoci={} candSeq={} acids={}) nucFrags(hetLoci={} candSeq={} nucs={})",
                refCoverageFragments.size(), refAminoAcidHetLoci.size(), candidateAminoAcidSequences.size(), refAminoAcids.size(),
                refNucleotideHetLoci.size(), candidateNucleotideSequences.size(), refNucleotides.size());

        List<FragmentAlleles> results = Lists.newArrayList();

        for(Fragment fragment : refCoverageFragments)
        {
            FragmentAlleles fragmentAlleles = create(
                    fragment, refAminoAcidHetLoci, candidateAminoAcidSequences, refAminoAcids,
                    refNucleotideHetLoci, candidateNucleotideSequences, refNucleotides);

            // drop wild-only alleles since their support can't be clearly established
            if(!fragmentAlleles.getFull().isEmpty())
                results.add(fragmentAlleles);
        }

        return results;
    }

    private static FragmentAlleles create(
            final Fragment fragment, final List<Integer> aminoAcidLoci, final List<HlaSequenceLoci> aminoAcidSequences,
            final List<Set<String>> refAminoAcids,
            final Map<String,List<Integer>> refNucleotideHetLoci, final List<HlaSequenceLoci> nucleotideSequences,
            final List<Set<String>> refNucleotides)
    {
        // look first for nucleotide support at the exon boundaries, then for amino acid support - full or wild
        Map<String,List<Integer>> fragNucleotideLociMap = Maps.newHashMap();

        refNucleotideHetLoci.entrySet().forEach(x -> fragNucleotideLociMap.put(
                x.getKey(), fragment.getNucleotideLoci().stream().filter(y -> x.getValue().contains(y)).collect(Collectors.toList())));

        final List<HlaAllele> fullNucleotideMatch = Lists.newArrayList();
        final List<HlaAllele> wildNucleotideMatch = Lists.newArrayList();

        Map<HlaAllele, SequenceMatchType> nucleotideAlleleMatches = findNucleotideMatches(
                fragment, refNucleotideHetLoci, nucleotideSequences, refNucleotides);

        fullNucleotideMatch.addAll(nucleotideAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList()));

        wildNucleotideMatch.addAll(nucleotideAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == WILD)
                .filter(x -> !fullNucleotideMatch.contains(x))
                .map(x -> x.getKey()).collect(Collectors.toList()));

        Map<HlaAllele, SequenceMatchType> aminoAcidAlleleMatches = findAminoAcidMatches(
                fragment, aminoAcidLoci, aminoAcidSequences, refAminoAcids);

        List<HlaAllele> fullAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> wildAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == WILD).map(x -> x.getKey()).collect(Collectors.toList());

        if(fullNucleotideMatch.isEmpty() && wildNucleotideMatch.isEmpty())
            return new FragmentAlleles(fragment, fullAminoAcidMatch, wildAminoAcidMatch);

        List<HlaAllele> consistentFull = fullAminoAcidMatch.stream()
                .filter(x -> fullNucleotideMatch.contains(x.asFourDigit())).collect(Collectors.toList());

        fullAminoAcidMatch.stream()
                .filter(x -> !wildAminoAcidMatch.contains(x))
                .filter(x -> wildNucleotideMatch.contains(x))
                .forEach(x -> wildAminoAcidMatch.add(x));

        // wildAminoAcidMatch.stream().filter(x -> !downgradedToWild.contains(x)).forEach(x -> do);

        return new FragmentAlleles(fragment, consistentFull, wildAminoAcidMatch);
    }

    private static Map<HlaAllele, SequenceMatchType> findNucleotideMatches(
            final Fragment fragment, final Map<String,List<Integer>> refNucleotideLociMap,
            final List<HlaSequenceLoci> nucleotideSequences, final List<Set<String>> refNucleotides)
    {
        Map<String,List<Integer>> fragNucleotideLociMap = Maps.newHashMap();

        Map<HlaAllele, SequenceMatchType> alleleMatches = Maps.newHashMap();

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

            // filter out wildcard bases at these exon boundaries
            List<Integer> alleleNucleotideLoci = filterWildcards(sequence, fragNucleotideLoci);

            SequenceMatchType matchType;
            if(alleleNucleotideLoci.isEmpty())
            {
                // if all bases being considered a wild, the treat it as full in nucleotide space and rely on the amino acid match type
                matchType = FULL;
            }
            else
            {
                String fragNucleotides = fragNucleotideSequences.get(allele.Gene);

                if(fragNucleotides == null)
                {
                    if(missedNucleotides.isEmpty())
                    {
                        fragNucleotides = fragment.nucleotides(alleleNucleotideLoci);
                    }
                    else
                    {
                        StringJoiner nucSj = new StringJoiner("");

                        for(Integer locus : alleleNucleotideLoci)
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

                matchType = sequence.determineMatchType(fragNucleotides, alleleNucleotideLoci);
            }

            if(matchType == SequenceMatchType.NONE)
                continue;

            // keep the best match
            Map.Entry<HlaAllele, SequenceMatchType> entryMatch = alleleMatches.entrySet().stream()
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

    private static Map<HlaAllele, SequenceMatchType> findAminoAcidMatches(
            final Fragment fragment, final List<Integer> aminoAcidLoci, final List<HlaSequenceLoci> aminoAcidSequences,
            final List<Set<String>> refAminoAcids)
    {
        Map<HlaAllele, SequenceMatchType> alleleMatches = Maps.newHashMap();

        List<Integer> fragAminoAcidLoci = fragment.getAminoAcidLoci().stream()
                .filter(x -> aminoAcidLoci.contains(x)).collect(Collectors.toList());

        if(fragAminoAcidLoci.isEmpty())
            return alleleMatches;

        // also attempt to retrive amino acids from low-qual nucleotides
        Map<Integer,String> missedAminoAcids = Maps.newHashMap();

        List<Integer> missedAminoAcidLoci = aminoAcidLoci.stream()
                .filter(x -> !fragAminoAcidLoci.contains(x)).collect(Collectors.toList());

        for(Integer missedLocus : missedAminoAcidLoci)
        {
            if(missedLocus >= refAminoAcids.size())
                continue;

            Set<String> candidateAminoAcids = refAminoAcids.get(missedLocus);

            String lowQualAminoAcid = fragment.getLowQualAminoAcid(missedLocus);

            if(!lowQualAminoAcid.isEmpty() && candidateAminoAcids.contains(lowQualAminoAcid))
            {
                fragAminoAcidLoci.add(missedLocus);
                missedAminoAcids.put(missedLocus, lowQualAminoAcid);
            }
        }

        Collections.sort(fragAminoAcidLoci);

        String fragmentAminoAcids;

        if(missedAminoAcids.isEmpty())
        {
            fragmentAminoAcids = fragment.aminoAcids(fragAminoAcidLoci);
        }
        else
        {
            StringJoiner aaSj = new StringJoiner("");

            for(Integer locus : fragAminoAcidLoci)
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
            HlaAllele allele = sequence.Allele;

            if(!fragment.getGenes().contains(allele.geneName()))
                continue;

            // ignore any wildcard loci at an exon boundary
            List<Integer> alleleFilteredLoci = filterExonBoundaryWildcards(sequence, fragAminoAcidLoci);

            String adjustedSequence = "";

            if(fragAminoAcidLoci.size() > alleleFilteredLoci.size())
            {
                for(int i = 0; i < fragAminoAcidLoci.size(); ++i)
                {
                    if(alleleFilteredLoci.contains(fragAminoAcidLoci.get(i)))
                        adjustedSequence += fragmentAminoAcids.charAt(i);
                }
            }
            else
            {
                adjustedSequence = fragmentAminoAcids;
            }

            SequenceMatchType matchType = sequence.determineMatchType(adjustedSequence, alleleFilteredLoci);
            if(matchType == SequenceMatchType.NONE)
                continue;

            Map.Entry<HlaAllele, SequenceMatchType> entryMatch = alleleMatches.entrySet().stream()
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

    public static void checkHlaYSupport(
            final String sampleId, final List<HlaSequenceLoci> hlaYSequences, List<Integer> refAminoAcidHetLoci,
            final List<FragmentAlleles> fragmentAlleles, final List<Fragment> refCoverageFragments)
    {
        // ignore fragments which don't contain any heterozygous locations
        int uniqueHlaY = 0;

        List<FragmentAlleles> matchedFragmentAlleles = Lists.newArrayList();

        for(Fragment fragment : refCoverageFragments)
        {
            List<Integer> fragAminoAcidLoci = fragment.getAminoAcidLoci().stream()
                    .filter(x -> refAminoAcidHetLoci.contains(x)).collect(Collectors.toList());

            if(fragAminoAcidLoci.isEmpty())
                continue;

            List<Integer> fragNucleotideLoci = fragment.getNucleotideLoci();

            boolean matchesY = false;
            FragmentAlleles matchedFrag = null;

            for(HlaSequenceLoci sequence : hlaYSequences)
            {
                String fragNucleotides = fragment.nucleotides(fragNucleotideLoci);

                SequenceMatchType matchType = sequence.determineMatchType(fragNucleotides, fragNucleotideLoci);
                if(matchType == SequenceMatchType.FULL)
                {
                    matchesY = true;

                    matchedFrag = fragmentAlleles.stream()
                            .filter(x -> x.getFragment().id().equals(fragment.id())).findFirst().orElse(null);

                    /*
                    if(LL_LOGGER.isDebugEnabled())
                    {
                        HlaAllele allele = sequence.Allele;

                        LL_LOGGER.debug("HLA-Y allele({}) fragment({}: {}) range({} -> {}) assignedGenes({})",
                                allele.toString(), fragment.id(), fragment.readInfo(),
                                fragment.getNucleotideLoci().get(0),
                                fragment.getNucleotideLoci().get(fragment.getNucleotideLoci().size() - 1),
                                matchedFrag != null ? matchedFrag.getFragment().getGenes() : "");
                    }
                    */

                    break;
                }
            }

            if(!matchesY)
                continue;

            if(matchedFrag == null)
            {
                ++uniqueHlaY;
            }
            else
            {
                matchedFragmentAlleles.add(matchedFrag);
            }
        }

        int totalHlaYFrags = uniqueHlaY  + matchedFragmentAlleles.size();
        double threshold = fragmentAlleles.size() * HLA_Y_FRAGMENT_THRESHOLD;
        boolean exceedsThreshold = uniqueHlaY >= threshold;

        if(totalHlaYFrags > 0)
        {
            LL_LOGGER.info("sample({}) HLA-Y fragments({} unique={}) shared={}) aboveThreshold({})",
                    sampleId, totalHlaYFrags, uniqueHlaY, matchedFragmentAlleles.size(), exceedsThreshold);

            if(exceedsThreshold)
            {
                matchedFragmentAlleles.forEach(x -> fragmentAlleles.remove(x));
            }
        }
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
                    .map(x -> new FragmentAlleles(x.getFragment(), Lists.newArrayList(stopLossAllele), Lists.newArrayList()))
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
