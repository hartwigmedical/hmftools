package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_Y_FRAGMENT_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_WILDCARD_FRAGMENTS;
import static com.hartwig.hmftools.lilac.LilacConstants.getAminoAcidExonBoundaries;
import static com.hartwig.hmftools.lilac.LilacConstants.getNucleotideExonBoundaries;
import static com.hartwig.hmftools.lilac.LilacConstants.longGeneName;
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
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

public class FragmentAlleleMapper
{
    private final Map<String,Map<Integer,Set<String>>> mGeneAminoAcidHetLociMap;
    private final Map<String,List<Integer>> mRefNucleotideHetLoci;
    private final List<Set<String>> mRefNucleotides;

    private int mStopLossFragments;
    private final List<HlaAllele> mStopLossAlleles;

    private final PerformanceCounter mPerfCounterFrag;

    public FragmentAlleleMapper(
            final Map<String, Map<Integer,Set<String>>> geneAminoAcidHetLociMap,
            final Map<String,List<Integer>> refNucleotideHetLoci, final List<Set<String>> refNucleotides)
    {
        mGeneAminoAcidHetLociMap = Maps.newHashMap();
        mGeneAminoAcidHetLociMap.putAll(geneAminoAcidHetLociMap);
        mRefNucleotideHetLoci = refNucleotideHetLoci;
        mRefNucleotides = refNucleotides;

        mStopLossFragments = 0;
        mStopLossAlleles = Lists.newArrayList();

        mPerfCounterFrag = new PerformanceCounter("Frags");
    }

    public void setHetAminoAcidLoci(final Map<String, Map<Integer,Set<String>>> geneAminoAcidHetLociMap)
    {
        mGeneAminoAcidHetLociMap.clear();
        mGeneAminoAcidHetLociMap.putAll(geneAminoAcidHetLociMap);
    }

    public void logPerfData()
    {
        mPerfCounterFrag.logStats();
    }

    public List<FragmentAlleles> createFragmentAlleles(
            final List<Fragment> refCoverageFragments, final List<HlaSequenceLoci> candidateAminoAcidSequences,
            final List<HlaSequenceLoci> candidateNucleotideSequences)
    {
        LL_LOGGER.info("building frag-alleles from aminoAcids(frags={} candSeq={}) nucFrags(hetLoci={} candSeq={} nucs={})",
                refCoverageFragments.size(), candidateAminoAcidSequences.size(),
                mRefNucleotideHetLoci.size(), candidateNucleotideSequences.size(), mRefNucleotides.size());

        mPerfCounterFrag.reset();

        List<FragmentAlleles> results = Lists.newArrayList();

        for(Fragment fragment : refCoverageFragments)
        {
            if(fragment.scope() == HLA_Y) // ignore if previously established
                continue;

            mPerfCounterFrag.start();
            FragmentAlleles fragAllele = mapFragmentToAlleles(fragment, candidateAminoAcidSequences, candidateNucleotideSequences);
            mPerfCounterFrag.stop();

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

        applyUniqueStopLossFragments(results);

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

            for(Map.Entry<String,Map<Integer,Set<String>>> geneEntry : mGeneAminoAcidHetLociMap.entrySet())
            {
                if(!fragment.getGenes().contains(longGeneName(geneEntry.getKey())))
                    continue;

                if(fragment.getAminoAcidLoci().stream().anyMatch(x -> geneEntry.getValue().containsKey(x)))
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
        Map<String,List<Integer>> fragNucleotideLociMap = Maps.newHashMap();

        mRefNucleotideHetLoci.entrySet().forEach(x -> fragNucleotideLociMap.put(
                x.getKey(), fragment.getNucleotideLoci().stream().filter(y -> x.getValue().contains(y)).collect(Collectors.toList())));

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
        List<HlaAllele> consistentFull = fullNucleotideMatch.stream()
                .filter(x -> fullAminoAcidMatch.contains(x) || homLociAminoAcidMatch.contains(x)).collect(Collectors.toList());

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
        Map<String,List<Integer>> fragGeneLociMap = Maps.newHashMap();
        Map<String,List<String>> fragGeneSequenceMap = Maps.newHashMap();

        Map<HlaAllele, SequenceMatchType> alleleMatches = Maps.newHashMap();

        // also attempt to retrieve amino acids from low-qual nucleotides
        Map<Integer,String> missedNucleotides = Maps.newHashMap();

        for(Map.Entry<String,List<Integer>> entry : mRefNucleotideHetLoci.entrySet())
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
                if(missedLocus >= mRefNucleotides.size())
                    continue;

                String lowQualNucleotide = fragment.getLowQualNucleotide(missedLocus);

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

            if(!fragment.getGenes().contains(allele.geneName()))
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

        Map<String,List<Integer>> fragGeneLociMap = Maps.newHashMap(); // per-gene map of heterozygous locations for this fragment
        Map<String,List<String>> fragGeneSequenceMap = Maps.newHashMap(); // per-gene map of fragment sequences at these het loci

        for(Map.Entry<String,Map<Integer,Set<String>>> geneEntry : mGeneAminoAcidHetLociMap.entrySet())
        {
            Map<Integer,Set<String>> hetLociSeqMap = geneEntry.getValue();

            List<Integer> fragAminoAcidLoci = fragment.getAminoAcidLoci().stream()
                    .filter(x -> hetLociSeqMap.containsKey(x)).collect(Collectors.toList());

            // also attempt to retrieve amino acids from low-qual nucleotides
            Map<Integer,String> missedAminoAcids = Maps.newHashMap();

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

            Collections.sort(fragAminoAcidLoci);

            final List<String> fragmentAminoAcids = Lists.newArrayListWithExpectedSize(fragAminoAcidLoci.size());

            for(Integer locus : fragAminoAcidLoci)
            {
                if(missedAminoAcids.containsKey(locus))
                    fragmentAminoAcids.add(missedAminoAcids.get(locus));
                else
                    fragmentAminoAcids.add(fragment.aminoAcid(locus));
            }

            fragGeneLociMap.put(geneEntry.getKey(), fragAminoAcidLoci);
            fragGeneSequenceMap.put(geneEntry.getKey(), fragmentAminoAcids);
        }

        for(HlaSequenceLoci sequence : aminoAcidSequences)
        {
            HlaAllele allele = sequence.Allele;

            if(!fragment.getGenes().contains(allele.geneName()))
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

    public boolean checkHlaYSupport(
            final List<HlaSequenceLoci> hlaYSequences, final List<FragmentAlleles> fragAlleles,
            final List<Fragment> fragments, boolean testThreshold)
    {
        // test for presence of HLA-Y and strip out from consideration any fragment mapping to it
        int uniqueHlaY = 0;

        List<FragmentAlleles> matchedFragmentAlleles = Lists.newArrayList();

        // ignore fragments which don't contain any heterozygous locations
        // only test heterozygous locations in A since HLA-Y matches its exon boundaries
        Set<Integer> aminoAcidHetLoci = mGeneAminoAcidHetLociMap.get(GENE_A).keySet();

        for(Fragment fragment : fragments)
        {
            List<Integer> fragAminoAcidLoci = fragment.getAminoAcidLoci().stream()
                    .filter(x -> aminoAcidHetLoci.contains(x)).collect(Collectors.toList());

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

                    matchedFrag = fragAlleles.stream()
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
                fragment.setScope(HLA_Y, true);
            }
            else
            {
                matchedFragmentAlleles.add(matchedFrag);
            }
        }

        int totalHlaYFrags = uniqueHlaY  + matchedFragmentAlleles.size();
        double threshold = fragments.size() * HLA_Y_FRAGMENT_THRESHOLD;
        boolean exceedsThreshold = uniqueHlaY >= threshold;

        if(totalHlaYFrags > 0)
        {
            if(testThreshold)
            {
                LL_LOGGER.info("HLA-Y fragments({} unique={}) shared={}) aboveThreshold({})",
                        totalHlaYFrags, uniqueHlaY, matchedFragmentAlleles.size(), exceedsThreshold);
            }

            if(exceedsThreshold || !testThreshold)
            {
                matchedFragmentAlleles.forEach(x -> fragAlleles.remove(x));
                matchedFragmentAlleles.forEach(x -> x.getFragment().setScope(HLA_Y, true));
            }
        }

        return exceedsThreshold;
    }

    public void setStopLossInfo(int stopLossFragments, final List<HlaAllele> stopLossAlleles)
    {
        mStopLossFragments = stopLossFragments;
        mStopLossAlleles.addAll(stopLossAlleles);
    }

    private void applyUniqueStopLossFragments(final List<FragmentAlleles> fragmentAlleles)
    {
        // create a unique (ie FULL-only) fragment allele for any confirmed known stop-loss allele
        if(mStopLossFragments == 0)
            return;

        for(HlaAllele stopLossAllele : mStopLossAlleles)
        {
            List<FragmentAlleles> sampleFragments = fragmentAlleles.stream()
                    .filter(x -> x.contains(stopLossAllele))
                    .map(x -> new FragmentAlleles(x.getFragment(), Lists.newArrayList(stopLossAllele), Lists.newArrayList()))
                    .collect(Collectors.toList());

            for(int i = 0; i < mStopLossFragments; ++i)
            {
                if(i >= sampleFragments.size())
                    break;

                fragmentAlleles.add(sampleFragments.get(i));
            }
        }
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

                LL_LOGGER.debug("wildcard allele({}) support from read({} {})",
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
        else
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
                fragAlleles.remove(index);
                ++removedFragAlleles;
            }
            else
            {
                ++index;
            }
        }

        LL_LOGGER.info("removed frags({}) with wildcard alleles removed", removedFragAlleles);
    }

    private static void removeUnsupportedWildcardAlleles(final List<HlaAllele> alleleList, final List<HlaAllele> unsupportedAlleles)
    {
        int index = 0;
        while(index < alleleList.size())
        {
            if(unsupportedAlleles.contains(alleleList.get(index)))
            {
                alleleList.remove(index);
            }
            else
            {
                ++index;
            }
        }
    }

}
