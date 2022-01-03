package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.LOG_UNMATCHED_HAPLOTYPE_SUPPORT;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import org.apache.commons.math3.util.Pair;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class HaplotypeQC
{
    public final int UnusedHaplotypes;
    public final int UnusedHaplotypeMaxFrags;
    public final int UnusedHaplotypesPon;

    public final List<Haplotype> UnmatchedHaplotypes;

    private static final List<Haplotype> PON_HAPLOTYPES = Lists.newArrayList();

    public HaplotypeQC(
            int unusedHaplotypes, int unusedHaplotypeMaxFrags, int unusedHaplotypesPon, final List<Haplotype> haplotypes)
    {
        UnusedHaplotypes = unusedHaplotypes;
        UnusedHaplotypeMaxFrags = unusedHaplotypeMaxFrags;
        UnusedHaplotypesPon = unusedHaplotypesPon;
        UnmatchedHaplotypes = haplotypes;
    }

    public List<String> header()
    {
        return Lists.newArrayList(
                "UnusedHaplotypes", "UnusedHaplotypeMaxFrags");
    }

    public List<String> body()
    {
        return Lists.newArrayList(
                String.valueOf(UnusedHaplotypes),
                String.valueOf(UnusedHaplotypeMaxFrags));
    }

    private static void loadPonHaplotypes()
    {
        if(!PON_HAPLOTYPES.isEmpty())
            return;

        final List<String> ponHaplotypes = new BufferedReader(new InputStreamReader(
                HaplotypeQC.class.getResourceAsStream("/pon/haplotypes.csv")))
                .lines().collect(Collectors.toList());

        ponHaplotypes.forEach(x -> PON_HAPLOTYPES.add(Haplotype.fromString(x)));
    }

    private static boolean inPon(final Haplotype haplotype)
    {
        return PON_HAPLOTYPES.stream().anyMatch(x -> x.contains(haplotype));
    }

    public static HaplotypeQC create(
            final List<HlaSequenceLoci> winners, final List<HlaSequenceLoci> hlaYSequences,
            final List<PhasedEvidence> evidence, final SequenceCount aminoAcidCount, final List<Fragment> unmatchedFragments)
    {
        loadPonHaplotypes();

        List<Haplotype> allUnmatched = Lists.newArrayList();
        evidence.stream()
                .map(x -> unmatchedHaplotype(x, LOG_UNMATCHED_HAPLOTYPE_SUPPORT, winners, aminoAcidCount, hlaYSequences))
                .forEach(x -> allUnmatched.addAll(x));

        Collections.sort(allUnmatched, new Haplotype.HaplotypeFragmentSorter());

        Set<String> haplotypeKeys = Sets.newHashSet();
        List<Haplotype> distinctHaplotypes = Lists.newArrayList();

        for(Haplotype haplotype : allUnmatched)
        {
            String key = String.format("%d-%d-%s", haplotype.StartLocus, haplotype.EndLocus, haplotype.Haplotype);
            if(haplotypeKeys.contains(key))
                continue;

            haplotypeKeys.add(key);
            distinctHaplotypes.add(haplotype);
        }

        findSupportingFragmentCounts(distinctHaplotypes, unmatchedFragments);

        Collections.sort(distinctHaplotypes, new Haplotype.HaplotypeStartLocusSorter());

        int pon = 0;
        int unusedCount = 0;
        int maxSupport = 0;

        for (Haplotype unmatched : distinctHaplotypes)
        {
            if(unmatched.matchingFragmentCount() < LOG_UNMATCHED_HAPLOTYPE_SUPPORT)
                continue;

            if (inPon(unmatched))
            {
                pon++;
                LL_LOGGER.info("  UNMATCHED_PON_HAPLTOYPE - {}", unmatched);
            }
            else
            {
                maxSupport = max(maxSupport, unmatched.matchingFragmentCount());
                unusedCount++;
                LL_LOGGER.warn("  UNMATCHED_HAPLTOYPE {}", unmatched);
            }
        }

        return new HaplotypeQC(unusedCount, maxSupport, pon, distinctHaplotypes);
    }

    private static boolean consistentWithAny(final PhasedEvidence phasedEvidence, final List<HlaSequenceLoci> winners, final String sequence)
    {
        return winners.stream().anyMatch(x -> x.consistentWith(sequence, phasedEvidence.getAminoAcidLoci()));
    }

    public static List<Haplotype> unmatchedHaplotype(
            final PhasedEvidence evidence, int minEvidence, final List<HlaSequenceLoci> winners, final SequenceCount aminoAcidCount,
            final List<HlaSequenceLoci> hlaYSequences)
    {
        // look through all phased evidence for AA sequences which are not supported by the winning alleles
        // ignore any wildcard sections
        Map<String,Integer> unmatched = evidence.getEvidence().entrySet().stream()
                .filter(x -> !consistentWithAny(evidence, winners, x.getKey()))
                .filter(x -> !consistentWithAny(evidence, hlaYSequences, x.getKey()))
                .filter(x -> x.getValue() >= minEvidence)
                .collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

        if(unmatched.isEmpty())
            return Lists.newArrayList();

        List<Haplotype> haplotypes = Lists.newArrayList();
        for(Map.Entry<String,Integer> entry : unmatched.entrySet())
        {
            Haplotype haplotype = Haplotype.create(
                    evidence.getAminoAcidLoci(), Pair.create(entry.getKey(), entry.getValue()), aminoAcidCount);

            haplotypes.add(haplotype);
        }

        return haplotypes;
    }

    private static void findSupportingFragmentCounts(final List<Haplotype> haplotypes, final List<Fragment> fragments)
    {
        for(Haplotype haplotype : haplotypes)
        {
            for(Fragment fragment : fragments)
            {
                if(!positionsOverlap(fragment.minAminoAcidLocus(), fragment.maxAminoAcidLocus(), haplotype.StartLocus, haplotype.EndLocus))
                    continue;

                boolean matches = true;
                int matchCount = 0;

                for(int locus = haplotype.StartLocus; locus <= haplotype.EndLocus; ++locus)
                {
                    int index = fragment.getAminoAcidLoci().indexOf(locus);
                    String fragmentAA = "";

                    if(index >= 0)
                    {
                        fragmentAA = fragment.getAminoAcids().get(index);
                    }
                    else
                    {
                        // retrieve the low-qual acid instead
                        fragmentAA = fragment.getLowQualAminoAcid(locus);

                        if(fragmentAA.isEmpty())
                            continue;
                    }

                    if(!fragmentAA.equals(haplotype.sequence(locus)))
                    {
                        matches = false;
                        break;
                    }

                    ++matchCount;
                }

                if(matches && matchCount > 0)
                {
                    //LL_LOGGER.debug("haplotype({}) supported by fragment({} {}) matchedAAs({})",
                    //        haplotype, fragment.id(), fragment.readInfo(), matchCount);

                    haplotype.addMatchingFragmentCount();
                }
            }
        }
    }
}
