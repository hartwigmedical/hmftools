package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_WILDCARD_FRAGMENTS;
import static com.hartwig.hmftools.lilac.ReferenceData.getAminoAcidExonBoundaries;
import static com.hartwig.hmftools.lilac.ReferenceData.getNucleotideExonBoundaries;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.NO_HET_LOCI;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.UNMATCHED_AMINO_ACID;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.WILD_ONLY;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.FULL;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.NO_LOCI;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.WILD;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.evidence.AminoAcid;
import com.hartwig.hmftools.lilac.evidence.Nucleotide;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

public class FragmentAlleleMapper
{
    private final Map<HlaGene, Map<Integer,Set<String>>> mGeneAminoAcidHetLociMap;
    private final Map<HlaGene, List<Integer>> mRefNucleotideHetLoci;
    private final List<Set<String>> mRefNucleotides;

    private final Map<HlaAllele,List<Fragment>> mStopLossAlleleFragments;

    public FragmentAlleleMapper(
            final Map<HlaGene, Map<Integer, Set<String>>> geneAminoAcidHetLociMap,
            final Map<HlaGene, List<Integer>> refNucleotideHetLoci, final List<Set<String>> refNucleotides)
    {
        mGeneAminoAcidHetLociMap = Maps.newHashMap();
        mGeneAminoAcidHetLociMap.putAll(geneAminoAcidHetLociMap);
        mRefNucleotideHetLoci = refNucleotideHetLoci;
        mRefNucleotides = refNucleotides;

        mStopLossAlleleFragments = Maps.newHashMap();
    }

    public void setHetAminoAcidLoci(final Map<HlaGene, Map<Integer, Set<String>>> geneAminoAcidHetLociMap)
    {
        mGeneAminoAcidHetLociMap.clear();
        mGeneAminoAcidHetLociMap.putAll(geneAminoAcidHetLociMap);
    }

    public void setKnownStopLossAlleleFragments(final Map<HlaAllele,List<Fragment>> knownStopLossFragments)
    {
        mStopLossAlleleFragments.putAll(knownStopLossFragments);
    }

    public List<FragmentAlleles> createFragmentAlleles(
            final List<Fragment> refCoverageFragments, final List<HlaSequenceLoci> candidateAminoAcidSequences,
            final List<HlaSequenceLoci> candidateNucleotideSequences)
    {
        LL_LOGGER.info("building frag-alleles from aminoAcids(frags={} candSeq={}) nucFrags(hetLoci={} candSeq={} nucs={}) knownIndels({})",
                refCoverageFragments.size(), candidateAminoAcidSequences.size(),
                mRefNucleotideHetLoci.size(), candidateNucleotideSequences.size(), mRefNucleotides.size(),
                mStopLossAlleleFragments.values().stream().mapToInt(x -> x.size()).sum());

        List<FragmentAlleles> results = Lists.newArrayList();

        for(Fragment fragment : refCoverageFragments)
        {
            if(fragment.scope() == HLA_Y) // ignore if previously established
                continue;

            FragmentAlleles fragAllele = checkStopLossAlleleFragments(fragment);

            if(fragAllele == null)
                fragAllele = mapFragmentToAlleles(fragment, candidateAminoAcidSequences, candidateNucleotideSequences);

            // drop wild-only alleles since their support can't be clearly established
            if(!fragAllele.getFull().isEmpty())
            {
                results.add(fragAllele);
            }
            else
            {
                setFragmentScope(fragAllele);
            }
        }

        return results;
    }

    private void setFragmentScope(final FragmentAlleles fragAllele)
    {
        Fragment fragment = fragAllele.getFragment();

        if(fragment.scope() == HLA_Y) // ignore if previously established
            return;

        fragment.clearScope();

        // will be set once solution is established
        if(!fragAllele.getFull().isEmpty())
            return;

        if(!fragAllele.getWild().isEmpty())
        {
            fragment.setScope(WILD_ONLY);
        }
        else
        {
            // was this fragment homozygous in the context of all genes
            boolean hasHetLoci = false;

            for(Map.Entry<HlaGene, Map<Integer, Set<String>>> geneEntry : mGeneAminoAcidHetLociMap.entrySet())
            {
                if(!fragment.genes().contains(geneEntry.getKey()))
                    continue;

                if(fragment.aminoAcidsByLoci()
                        .values()
                        .stream()
                        .mapToInt(AminoAcid::locus)
                        .anyMatch(x -> geneEntry.getValue().containsKey(x)))
                {
                    hasHetLoci = true;
                    break;
                }
            }

            if(!hasHetLoci)
            {
                fragment.setScope(NO_HET_LOCI);
            }
            else
            {
                fragment.setScope(UNMATCHED_AMINO_ACID);
            }
        }
    }

    private FragmentAlleles mapFragmentToAlleles(
            final Fragment fragment, final List<HlaSequenceLoci> aminoAcidSequences, final List<HlaSequenceLoci> nucleotideSequences)
    {
        // look first for nucleotide support at the exon boundaries, then for amino acid support - full or wild
        Map<HlaGene, List<Integer>> fragNucleotideLociMap = Maps.newHashMap();

        mRefNucleotideHetLoci.entrySet().forEach(x -> fragNucleotideLociMap.put(x.getKey(), fragment.nucleotidesByLoci()
                        .values()
                        .stream()
                        .map(Nucleotide::locus)
                        .filter(y -> x.getValue().contains(y))
                        .collect(Collectors.toList())));

        Map<HlaAllele, SequenceMatchType> nucleotideAlleleMatches = findNucleotideMatches(fragment, nucleotideSequences);

        List<HlaAllele> fullNucleotideMatch = nucleotideAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> wildNucleotideMatch = nucleotideAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == WILD)
                .map(x -> x.getKey()).collect(Collectors.toList());

        Map<HlaAllele, SequenceMatchType> aminoAcidAlleleMatches = findAminoAcidMatches(fragment, aminoAcidSequences);

        List<HlaAllele> fullAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> wildAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == WILD).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> homLociAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == NO_LOCI).map(x -> x.getKey()).collect(Collectors.toList());

        if(fullNucleotideMatch.isEmpty() && wildNucleotideMatch.isEmpty())
        {
            // do not allow wild-only (ie no full) if there are homozygous matches
            if(fullAminoAcidMatch.isEmpty() && !wildAminoAcidMatch.isEmpty() && !homLociAminoAcidMatch.isEmpty())
                return new FragmentAlleles(fragment, Lists.newArrayList(), Lists.newArrayList());

            return new FragmentAlleles(fragment, fullAminoAcidMatch, wildAminoAcidMatch);
        }

        // otherwise look for matching nuc and amino-acid full matches
        Set<HlaAllele> fullAminoAcidMatchSet = Sets.newHashSet(fullAminoAcidMatch);
        Set<HlaAllele> homLociAminoAcidMatchSet = Sets.newHashSet(homLociAminoAcidMatch);
        List<HlaAllele> consistentFull = fullNucleotideMatch.stream()
                .filter(x -> fullAminoAcidMatchSet.contains(x) || homLociAminoAcidMatchSet.contains(x)).collect(Collectors.toList());

        // otherwise down-grade the full matches to wild
        fullAminoAcidMatch.stream()
                .filter(x -> !wildAminoAcidMatch.contains(x))
                .filter(x -> wildNucleotideMatch.contains(x))
                .forEach(x -> wildAminoAcidMatch.add(x));

        return new FragmentAlleles(fragment, consistentFull, wildAminoAcidMatch);
    }

    private Map<HlaAllele, SequenceMatchType> findNucleotideMatches(
            final Fragment fragment, final List<HlaSequenceLoci> nucleotideSequences)
    {
        Map<HlaGene, List<Integer>> fragGeneLociMap = Maps.newHashMap();
        Map<HlaGene, List<String>> fragGeneSequenceMap = Maps.newHashMap();

        Map<HlaAllele, SequenceMatchType> alleleMatches = Maps.newHashMap();

        // also attempt to retrieve amino acids from low-qual nucleotides
        Map<Integer, String> missedNucleotides = Maps.newHashMap();

        for(Map.Entry<HlaGene, List<Integer>> entry : mRefNucleotideHetLoci.entrySet())
        {
            HlaGene gene = entry.getKey();
            List<Integer> refNucleotideLoci = entry.getValue();

            List<Integer> fragmentMatchedLoci = fragment.nucleotidesByLoci()
                    .values()
                    .stream()
                    .map(Nucleotide::locus)
                    .filter(y -> refNucleotideLoci.contains(y))
                    .collect(Collectors.toList());

            // also check for support from low-qual reads at this same location
            List<Integer> missedNucleotideLoci = refNucleotideLoci.stream()
                    .filter(x -> !fragment.nucleotidesByLoci().containsKey(x)).collect(Collectors.toList());

            for(Integer missedLocus : missedNucleotideLoci)
            {
                if(missedLocus >= mRefNucleotides.size())
                    continue;

                String lowQualNucleotide = fragment.getRawNucleotide(missedLocus);

                if(lowQualNucleotide.isEmpty())
                    continue;

                Set<String> candidateNucleotides = mRefNucleotides.get(missedLocus);

                if(candidateNucleotides.contains(lowQualNucleotide))
                {
                    fragmentMatchedLoci.add(missedLocus);
                    missedNucleotides.put(missedLocus, lowQualNucleotide);
                }
            }

            Collections.sort(fragmentMatchedLoci);

            final List<String> fragmentNucleotides = Lists.newArrayListWithExpectedSize(fragmentMatchedLoci.size());

            for(Integer locus : fragmentMatchedLoci)
            {
                if(missedNucleotides.containsKey(locus))
                    fragmentNucleotides.add(missedNucleotides.get(locus));
                else
                    fragmentNucleotides.add(fragment.nucleotide(locus));
            }

            fragGeneLociMap.put(gene, fragmentMatchedLoci);
            fragGeneSequenceMap.put(gene, fragmentNucleotides);
        }

        if(fragGeneLociMap.values().stream().allMatch(x -> x.isEmpty()))
            return alleleMatches;

        for(HlaSequenceLoci sequence : nucleotideSequences)
        {
            HlaAllele allele = sequence.Allele;

            HlaAllele proteinAllele = allele.asFourDigit();

            SequenceMatchType existingMatch = alleleMatches.get(proteinAllele);

            if(existingMatch != null && existingMatch == FULL)
                continue;

            if(!fragment.genes().contains(allele.Gene))
                continue;

            List<Integer> fragNucleotideLoci = fragGeneLociMap.get(allele.Gene);
            if(fragNucleotideLoci.isEmpty())
                continue;

            // filter out wildcard bases at these exon boundaries
            List<String> fragmentNucleotides = fragGeneSequenceMap.get(allele.Gene);

            if(sequence.hasExonBoundaryWildcards())
            {
                // mustn't change the original
                fragmentNucleotides = fragmentNucleotides.stream().collect(Collectors.toList());
                fragNucleotideLoci = fragNucleotideLoci.stream().collect(Collectors.toList());

                // ignore any wildcard loci at an exon boundary
                List<Integer> nucleotideExonBoundaries = getNucleotideExonBoundaries(sequence.Allele.Gene);

                int index = 0;
                while(index < fragNucleotideLoci.size())
                {
                    int locus = fragNucleotideLoci.get(index);
                    boolean wildcardExonBoundary = nucleotideExonBoundaries.contains(locus)
                            && locus < sequence.length() && sequence.sequence(locus).equals(WILD_STR);

                    if(!wildcardExonBoundary)
                    {
                        ++index;
                    }
                    else
                    {
                        fragNucleotideLoci.remove(index);
                        fragmentNucleotides.remove(index);
                    }
                }
            }

            SequenceMatchType matchType;
            if(fragNucleotideLoci.isEmpty())
            {
                // if all bases being considered a wild, the treat it as full in nucleotide space and rely on the amino acid match type
                matchType = FULL;
            }
            else
            {
                matchType = sequence.determineMatchType(fragmentNucleotides, fragNucleotideLoci);
            }

            if(matchType == SequenceMatchType.MISMATCH)
                continue;

            if(existingMatch == null || matchType.isBetter(existingMatch))
            {
                alleleMatches.put(proteinAllele, matchType);
            }
        }

        return alleleMatches;
    }

    private Map<HlaAllele, SequenceMatchType> findAminoAcidMatches(final Fragment fragment, final List<HlaSequenceLoci> aminoAcidSequences)
    {
        Map<HlaAllele,SequenceMatchType> alleleMatches = Maps.newHashMap();

        Map<HlaGene, List<Integer>> fragGeneLociMap = Maps.newHashMap(); // per-gene map of heterozygous locations for this fragment
        Map<HlaGene, List<String>> fragGeneSequenceMap = Maps.newHashMap(); // per-gene map of fragment sequences at these het loci

        for(Map.Entry<HlaGene, Map<Integer, Set<String>>> geneEntry : mGeneAminoAcidHetLociMap.entrySet())
        {
            Map<Integer, Set<String>> hetLociSeqMap = geneEntry.getValue();

            NavigableSet<Integer> fragAminoAcidLoci = fragment.aminoAcidsByLoci().keySet().stream()
                    .filter(x -> hetLociSeqMap.containsKey(x))
                    .collect(Collectors.toCollection(Sets::newTreeSet));

            // also attempt to retrieve amino acids from low-qual nucleotides
            Map<Integer, String> missedAminoAcids = Maps.newHashMap();

            List<Integer> missedAminoAcidLoci = hetLociSeqMap.keySet().stream()
                    .filter(x -> !fragAminoAcidLoci.contains(x)).collect(Collectors.toList());

            for(Integer missedLocus : missedAminoAcidLoci)
            {
                String lowQualAminoAcid = fragment.getLowQualAminoAcid(missedLocus);

                if(lowQualAminoAcid.isEmpty())
                    continue;

                Set<String> candidateAminoAcids = hetLociSeqMap.get(missedLocus);

                if(candidateAminoAcids.contains(lowQualAminoAcid))
                {
                    fragAminoAcidLoci.add(missedLocus);
                    missedAminoAcids.put(missedLocus, lowQualAminoAcid);
                }
            }

            final List<String> fragmentAminoAcids = Lists.newArrayListWithExpectedSize(fragAminoAcidLoci.size());

            for(Integer locus : fragAminoAcidLoci)
            {
                if(missedAminoAcids.containsKey(locus))
                    fragmentAminoAcids.add(missedAminoAcids.get(locus));
                else
                    fragmentAminoAcids.add(fragment.aminoAcid(locus));
            }

            fragGeneLociMap.put(geneEntry.getKey(), Lists.newArrayList(fragAminoAcidLoci));
            fragGeneSequenceMap.put(geneEntry.getKey(), fragmentAminoAcids);
        }

        for(HlaSequenceLoci sequence : aminoAcidSequences)
        {
            HlaAllele allele = sequence.Allele;

            if(!fragment.genes().contains(allele.Gene))
            {
                alleleMatches.put(allele, NO_LOCI);
                continue;
            }

            List<Integer> fragAminoAcidLoci = fragGeneLociMap.get(allele.Gene);

            if(fragAminoAcidLoci == null || fragAminoAcidLoci.isEmpty()) // not supported by this gene or homozygous in all locations
            {
                alleleMatches.put(allele, NO_LOCI);
                continue;
            }

            List<String> fragmentAminoAcids = fragGeneSequenceMap.get(allele.Gene);

            if(sequence.hasExonBoundaryWildcards())
            {
                // mustn't change the original
                fragmentAminoAcids = fragmentAminoAcids.stream().collect(Collectors.toList());
                fragAminoAcidLoci = fragAminoAcidLoci.stream().collect(Collectors.toList());

                // ignore any wildcard loci at an exon boundary
                List<Integer> aminoAcidExonBoundaries = getAminoAcidExonBoundaries(sequence.Allele.Gene);

                int index = 0;
                while(index < fragAminoAcidLoci.size())
                {
                    int locus = fragAminoAcidLoci.get(index);
                    boolean wildcardExonBoundary = aminoAcidExonBoundaries.contains(locus)
                            && locus < sequence.length() && sequence.sequence(locus).equals(WILD_STR);

                    if(!wildcardExonBoundary)
                    {
                        ++index;
                    }
                    else
                    {
                        fragAminoAcidLoci.remove(index);
                        fragmentAminoAcids.remove(index);
                    }
                }
            }

            SequenceMatchType matchType = sequence.determineMatchType(fragmentAminoAcids, fragAminoAcidLoci);
            if(matchType == SequenceMatchType.MISMATCH)
                continue;

            alleleMatches.put(allele, matchType);
        }

        return alleleMatches;
    }

    private FragmentAlleles checkStopLossAlleleFragments(Fragment fragment)
    {
        for(Map.Entry<HlaAllele,List<Fragment>> entry : mStopLossAlleleFragments.entrySet())
        {
            if(entry.getValue().contains(fragment))
            {
                return new FragmentAlleles(fragment, Lists.newArrayList(entry.getKey()), Lists.newArrayList());
            }
        }

        return null;
    }

    public static Set<HlaAllele> findWildcardAlleles(final List<FragmentAlleles> fragAlleles)
    {
        Set<HlaAllele> wildcardAlleles = Sets.newHashSet();

        for(FragmentAlleles fragAllele : fragAlleles)
        {
            fragAllele.getFull().stream().filter(x -> x.hasWildcards()).forEach(x -> wildcardAlleles.add(x));
            fragAllele.getWild().stream().filter(x -> x.hasWildcards()).forEach(x -> wildcardAlleles.add(x));
        }

        return wildcardAlleles;
    }

    public static List<HlaAllele> findUnsupportedWildcards(final List<FragmentAlleles> fragAlleles, final Set<HlaAllele> wildcardAlleles)
    {
        // across all the candidates, find wildcards which have 2+ non-wild fragments which don't support a non-wildcard allele
        // remove unsupported wildcard alleles from fragments and return the list of supported wildcard alleles
        List<HlaAllele> supportedAlleles = Lists.newArrayList();

        for(HlaAllele allele : wildcardAlleles)
        {
            int supportCount = 0;

            for(FragmentAlleles fragAllele : fragAlleles)
            {
                if(!fragAllele.getFull().contains(allele))
                    continue;

                // ignore fragments supporting a non-wildcard allele
                if(fragAllele.getFull().stream().anyMatch(x -> !x.hasWildcards()))
                    continue;

                ++supportCount;

                LL_LOGGER.trace("wildcard allele({}) support from read({} {})",
                        allele, fragAllele.getFragment().id(), fragAllele.getFragment().readInfo());
            }

            if(supportCount >= MIN_WILDCARD_FRAGMENTS)
            {
                supportedAlleles.add(allele);
            }
        }

        if(!supportedAlleles.isEmpty())
        {
            LL_LOGGER.info("filtered wildcards({} -> {}) supported alleles: {}",
                    wildcardAlleles.size(), supportedAlleles.size(), HlaAllele.toString(supportedAlleles));
        }
        else if(!wildcardAlleles.isEmpty())
        {
            LL_LOGGER.info("removed all {} wildcard candidates", wildcardAlleles.size());
        }

        return wildcardAlleles.stream().filter(x -> !supportedAlleles.contains(x)).collect(Collectors.toList());
    }

    public static void filterUnsupportedWildcardFragments(final List<FragmentAlleles> fragAlleles, final List<HlaAllele> unsupportedAlleles)
    {
        // remove these alleles from fragment alleles
        int index = 0;
        int removedFragAlleles = 0;

        while(index < fragAlleles.size())
        {
            FragmentAlleles fragAllele = fragAlleles.get(index);

            removeUnsupportedWildcardAlleles(fragAllele.getFull(), unsupportedAlleles);
            removeUnsupportedWildcardAlleles(fragAllele.getWild(), unsupportedAlleles);

            if(fragAllele.getFull().isEmpty() && fragAllele.getWild().isEmpty())
            {
                fragAllele.getFragment().setScope(WILD_ONLY);
                fragAlleles.remove(index);
                ++removedFragAlleles;
            }
            else
            {
                ++index;
            }
        }

        if(removedFragAlleles > 0)
        {
            LL_LOGGER.debug("removed frags({}) with wildcard alleles removed", removedFragAlleles);
        }
    }

    private static void removeUnsupportedWildcardAlleles(final Set<HlaAllele> alleleList, final List<HlaAllele> unsupportedAlleles)
    {
        alleleList.removeAll(unsupportedAlleles);
    }
}
