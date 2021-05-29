package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_Y_FRAGMENT_THRESHOLD;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.CANDIDATE;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.HLA_Y;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.NO_HET_LOCI;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.UNMATCHED_AMINO_ACID;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.filterExonBoundaryWildcards;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceLoci.filterWildcards;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.FULL;
import static com.hartwig.hmftools.lilac.seq.SequenceMatchType.WILD;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceMatchType;

public class FragmentAlleleMapper
{
    private final Map<String, Map<Integer,List<String>>> mGeneAminoAcidHetLociMap;
    private final Map<String,List<Integer>> mRefNucleotideHetLoci;
    private final List<Set<String>> mRefNucleotides;

    private final PerformanceCounter mPerfCounterFrag;
    private final PerformanceCounter mPerfCounterAcid;
    private final PerformanceCounter mPerfCounterNuc;

    public FragmentAlleleMapper(
            final Map<String, Map<Integer,List<String>>> geneAminoAcidHetLociMap,
            final Map<String,List<Integer>> refNucleotideHetLoci, final List<Set<String>> refNucleotides)
    {
        mGeneAminoAcidHetLociMap = geneAminoAcidHetLociMap;
        mRefNucleotideHetLoci = refNucleotideHetLoci;
        mRefNucleotides = refNucleotides;

        mPerfCounterFrag = new PerformanceCounter("Frags");
        mPerfCounterNuc = new PerformanceCounter("Nuc");
        mPerfCounterAcid = new PerformanceCounter("AminoAcid");
    }

    public void logPerfData()
    {
        mPerfCounterFrag.logStats();
        mPerfCounterNuc.logStats();
        mPerfCounterAcid.logStats();
    }

    public List<FragmentAlleles> createFragmentAlleles(
            final List<Fragment> refCoverageFragments, final List<HlaSequenceLoci> candidateAminoAcidSequences,
            final List<HlaSequenceLoci> candidateNucleotideSequences)
    {
        LL_LOGGER.info("building frag-alleles from aminoAcids(frags={} candSeq={}) nucFrags(hetLoci={} candSeq={} nucs={})",
                refCoverageFragments.size(), candidateAminoAcidSequences.size(),
                mRefNucleotideHetLoci.size(), candidateNucleotideSequences.size(), mRefNucleotides.size());

        List<FragmentAlleles> results = Lists.newArrayList();

        for(Fragment fragment : refCoverageFragments)
        {
            mPerfCounterFrag.start();
            FragmentAlleles fragmentAlleles = mapFragmentToAlleles(fragment, candidateAminoAcidSequences, candidateNucleotideSequences);
            mPerfCounterFrag.stop();

            // drop wild-only alleles since their support can't be clearly established
            if(!fragmentAlleles.getFull().isEmpty())
            {
                results.add(fragmentAlleles);
            }
            else
            {
                // LL_LOGGER.debug("frag({}: {}) unassigned", fragment.id(), fragment.readInfo());

                if(!fragmentAlleles.getWild().isEmpty())
                {
                    fragment.setScope(CANDIDATE);
                }
                else
                {
                    if(fragment.getAminoAcidLoci().stream()
                            .noneMatch(x -> mGeneAminoAcidHetLociMap.values().stream().anyMatch(y -> y.containsKey(x))))
                    {
                        fragment.setScope(NO_HET_LOCI);
                    }
                    else
                    {
                        fragment.setScope(UNMATCHED_AMINO_ACID);
                    }
                }
            }
        }

        return results;
    }

    private FragmentAlleles mapFragmentToAlleles(
            final Fragment fragment, final List<HlaSequenceLoci> aminoAcidSequences, final List<HlaSequenceLoci> nucleotideSequences)
    {
        // look first for nucleotide support at the exon boundaries, then for amino acid support - full or wild
        Map<String,List<Integer>> fragNucleotideLociMap = Maps.newHashMap();

        mRefNucleotideHetLoci.entrySet().forEach(x -> fragNucleotideLociMap.put(
                x.getKey(), fragment.getNucleotideLoci().stream().filter(y -> x.getValue().contains(y)).collect(Collectors.toList())));

        mPerfCounterNuc.start();
        Map<HlaAllele, SequenceMatchType> nucleotideAlleleMatches = findNucleotideMatches(fragment, nucleotideSequences);
        mPerfCounterNuc.stop();

        List<HlaAllele> fullNucleotideMatch = nucleotideAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> wildNucleotideMatch = nucleotideAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == WILD)
                .map(x -> x.getKey()).collect(Collectors.toList());

        mPerfCounterAcid.start();
        Map<HlaAllele, SequenceMatchType> aminoAcidAlleleMatches = findAminoAcidMatches(fragment, aminoAcidSequences);
        mPerfCounterAcid.stop();

        List<HlaAllele> fullAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == FULL).map(x -> x.getKey()).collect(Collectors.toList());

        List<HlaAllele> wildAminoAcidMatch = aminoAcidAlleleMatches.entrySet().stream()
                .filter(x -> x.getValue() == WILD).map(x -> x.getKey()).collect(Collectors.toList());

        if(fullNucleotideMatch.isEmpty() && wildNucleotideMatch.isEmpty())
            return new FragmentAlleles(fragment, fullAminoAcidMatch, wildAminoAcidMatch);

        List<HlaAllele> consistentFull = fullAminoAcidMatch.stream()
                .filter(x -> fullNucleotideMatch.contains(x)).collect(Collectors.toList());

        fullAminoAcidMatch.stream()
                .filter(x -> !wildAminoAcidMatch.contains(x))
                .filter(x -> wildNucleotideMatch.contains(x))
                .forEach(x -> wildAminoAcidMatch.add(x));

        return new FragmentAlleles(fragment, consistentFull, wildAminoAcidMatch);
    }

    private Map<HlaAllele, SequenceMatchType> findNucleotideMatches(
            final Fragment fragment, final List<HlaSequenceLoci> nucleotideSequences)
    {
        Map<String,List<Integer>> fragNucleotideLociMap = Maps.newHashMap();

        Map<HlaAllele, SequenceMatchType> alleleMatches = Maps.newHashMap();

        // also attempt to retrive amino acids from low-qual nucleotides
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

            fragNucleotideLociMap.put(gene, fragmentMatchedLoci);
        }

        if(fragNucleotideLociMap.values().stream().allMatch(x -> x.isEmpty()))
            return alleleMatches;

        fragNucleotideLociMap.values().forEach(x -> Collections.sort(x));

        Map<String,String> fragNucleotideSequences = Maps.newHashMap();

        for(HlaSequenceLoci sequence : nucleotideSequences)
        {
            HlaAllele allele = sequence.Allele;

            HlaAllele proteinAllele = allele.asFourDigit();

            SequenceMatchType existingMatch = alleleMatches.get(proteinAllele);

            if(existingMatch != null && existingMatch == FULL)
                continue;

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

            if(existingMatch == null || matchType.isBetter(existingMatch))
            {
                alleleMatches.put(proteinAllele, matchType);
            }

            /*
            // keep or improve the best match
            Map.Entry<HlaAllele, SequenceMatchType> entryMatch = alleleMatches.entrySet().stream()
                    .filter(x -> x.getKey() == allele.asFourDigit()).findFirst().orElse(null);

            if(entryMatch == null || matchType.isBetter(entryMatch.getValue()))
            {
                if(entryMatch != null)
                    alleleMatches.remove(entryMatch.getKey());

                alleleMatches.put(allele.asFourDigit(), matchType);
            }
            */
        }

        return alleleMatches;
    }

    private Map<HlaAllele, SequenceMatchType> findAminoAcidMatches(final Fragment fragment, final List<HlaSequenceLoci> aminoAcidSequences)
    {
        Map<HlaAllele,SequenceMatchType> alleleMatches = Maps.newHashMap();

        Map<String,List<Integer>> fragGeneLociMap = Maps.newHashMap(); // per-gene map of heterozygous locations for this fragment
        Map<String,String> fragGeneSequenceMap = Maps.newHashMap(); // per-gene map of fragment sequences at these het loci

        for(Map.Entry<String,Map<Integer,List<String>>> geneEntry : mGeneAminoAcidHetLociMap.entrySet())
        {
            Map<Integer,List<String>> hetLociSeqMap = geneEntry.getValue();

            List<Integer> fragAminoAcidLoci = fragment.getAminoAcidLoci().stream()
                    .filter(x -> hetLociSeqMap.containsKey(x)).collect(Collectors.toList());

            if(fragAminoAcidLoci.isEmpty())
                continue;

            // also attempt to retrieve amino acids from low-qual nucleotides
            Map<Integer, String> missedAminoAcids = Maps.newHashMap();

            List<Integer> missedAminoAcidLoci = hetLociSeqMap.keySet().stream()
                    .filter(x -> !fragAminoAcidLoci.contains(x)).collect(Collectors.toList());

            for(Integer missedLocus : missedAminoAcidLoci)
            {
                String lowQualAminoAcid = fragment.getLowQualAminoAcid(missedLocus);

                if(lowQualAminoAcid.isEmpty())
                    continue;

                List<String> candidateAminoAcids = hetLociSeqMap.get(missedLocus);

                if(candidateAminoAcids.contains(lowQualAminoAcid))
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

            fragGeneLociMap.put(geneEntry.getKey(), fragAminoAcidLoci);
            fragGeneSequenceMap.put(geneEntry.getKey(), fragmentAminoAcids);
        }

        for(HlaSequenceLoci sequence : aminoAcidSequences)
        {
            HlaAllele allele = sequence.Allele;

            if(!fragment.getGenes().contains(allele.geneName()))
                continue;

            List<Integer> fragAminoAcidLoci = fragGeneLociMap.get(allele.Gene);

            if(fragAminoAcidLoci == null) // not supported by this gene or homozygous in all locations
                continue;

            String fragmentAminoAcids = fragGeneSequenceMap.get(allele.Gene);

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

            alleleMatches.put(allele, matchType);

            /*
            Map.Entry<HlaAllele, SequenceMatchType> entryMatch = alleleMatches.entrySet().stream()
                    .filter(x -> x.getKey() == allele).findFirst().orElse(null);

            if(entryMatch == null || matchType.isBetter(entryMatch.getValue()))
            {
                if(entryMatch != null)
                    alleleMatches.remove(entryMatch.getKey());
            }
             */
        }

        return alleleMatches;
    }

    public boolean checkHlaYSupport(
            final List<HlaSequenceLoci> hlaYSequences, final List<FragmentAlleles> fragAlleles, final List<Fragment> fragments)
    {
        // ignore fragments which don't contain any heterozygous locations
        int uniqueHlaY = 0;

        List<FragmentAlleles> matchedFragmentAlleles = Lists.newArrayList();

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
                fragment.setScope(HLA_Y);
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
            LL_LOGGER.info("HLA-Y fragments({} unique={}) shared={}) aboveThreshold({})",
                    totalHlaYFrags, uniqueHlaY, matchedFragmentAlleles.size(), exceedsThreshold);

            if(exceedsThreshold)
            {
                matchedFragmentAlleles.forEach(x -> fragAlleles.remove(x));
                matchedFragmentAlleles.forEach(x -> x.getFragment().setScope(HLA_Y, true));
            }
        }

        return exceedsThreshold;
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
