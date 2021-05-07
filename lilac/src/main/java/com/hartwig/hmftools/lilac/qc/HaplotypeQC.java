package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.max;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.sun.tools.javac.util.Pair;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class HaplotypeQC
{
    private final int mUnusedHaplotypes;
    private final int mUnusedHaplotypeMaxSupport;
    private final int mUnusedHaplotypeMaxLength;
    private final int mUnusedHaplotypesPon;

    private static final int MIN_SUPPORT = 7;

    private static final List<Haplotype> PON_HAPLOTYPES = Lists.newArrayList();

    public final List<String> header()
    {
        return Lists.newArrayList(
                "unusedHaplotypes", "unusedHaplotypeMaxSupport", "unusedHaplotypeMaxLength", "unusedHaplotypesPon");
    }

    public final List<String> body()
    {
        return Lists.newArrayList(
                String.valueOf(mUnusedHaplotypes),
                String.valueOf(mUnusedHaplotypeMaxSupport), 
                String.valueOf(mUnusedHaplotypeMaxLength),
                String.valueOf(mUnusedHaplotypesPon));
    }

    public final int getUnusedHaplotypes()
    {
        return mUnusedHaplotypes;
    }

    public HaplotypeQC(int unusedHaplotypes, int unusedHaplotypeMaxSupport, int unusedHaplotypeMaxLength, int unusedHaplotypesPon)
    {
        mUnusedHaplotypes = unusedHaplotypes;
        mUnusedHaplotypeMaxSupport = unusedHaplotypeMaxSupport;
        mUnusedHaplotypeMaxLength = unusedHaplotypeMaxLength;
        mUnusedHaplotypesPon = unusedHaplotypesPon;

        final List<String> ponHaplotypes = new BufferedReader(new InputStreamReader(
                HaplotypeQC.class.getResourceAsStream("/pon/haplotypes.txt")))
                .lines().collect(Collectors.toList());

        ponHaplotypes.forEach(x -> PON_HAPLOTYPES.add(Haplotype.fromString(x)));
    }

    private static boolean inPon(final Haplotype haplotype)
    {
        return PON_HAPLOTYPES.contains(haplotype);
    }

    public static HaplotypeQC create(
            int minEvidence, final List<HlaSequenceLoci> winners, final List<PhasedEvidence> evidence, final SequenceCount aminoAcidCount)
    {
        List<Haplotype> allUnmatched = Lists.newArrayList();
        evidence.stream().map(x -> unmatchedHaplotype(x, minEvidence, winners, aminoAcidCount)).forEach(x -> allUnmatched.addAll(x));
        Collections.sort(allUnmatched, new Haplotype.HaplotypeFragmentSorter());

        Map<String,Haplotype> distinctHaplotypeMap = allUnmatched.stream()
                .collect(Collectors.toMap(entry -> String.format("%d-%d-%s", entry.StartLocus, entry.EndLocus, entry.Haplotype), entry -> entry));

        List<Haplotype> distinctHaplotypes = distinctHaplotypeMap.values().stream().collect(Collectors.toList());
        Collections.sort(distinctHaplotypes, new Haplotype.HaplotypeStartLocusSorter());

        /*
        val allUnmatched = evidence
                .flatMap { it.unmatchedHaplotype(minEvidence, winners, aminoAcidCount) }
                    .sortedByDescending { it.supportingFragments }
                    .distinctBy { x -> (x.startLocus.toString() + x.endLocus.toString() + x.haplotype) }
                    .sortedBy { it.startLocus }
         */

        int pon = 0;
        int unusedCount = 0;
        int maxSupport = 0;
        int maxLength = 0;

        for (Haplotype unmatched : allUnmatched)
        {
            if (unmatched.SupportingFragments >= MIN_SUPPORT)
            {
                if (inPon(unmatched))
                {
                    pon++;
                    LL_LOGGER.info("  UNMATCHED_PON_HAPLTOYPE - {}", unmatched);
                }
                else
                {
                    maxSupport = max(maxSupport, unmatched.SupportingFragments);
                    maxLength = max(maxLength, unmatched.Haplotype.length());
                    unusedCount++;

                    LL_LOGGER.warn("  UNMATCHED_HAPLTOYPE - $unmatched");
                }
            }
        }

        return new HaplotypeQC(unusedCount, maxSupport, maxLength, pon);
    }

    private static boolean consistentWithAny(final PhasedEvidence phasedEvidence, final List<HlaSequenceLoci> winners, final String sequence)
    {
        return winners.stream().anyMatch(x -> x.consistentWith(sequence, phasedEvidence.getAminoAcidIndices()));
    }

    public static List<Haplotype> unmatchedHaplotype(
            final PhasedEvidence evidence, int minEvidence,
            final List<HlaSequenceLoci> winners, final SequenceCount aminoAcidCount)
    {
        Map<String,Integer> unmatched = evidence.getEvidence().entrySet().stream()
                .filter(x -> !consistentWithAny(evidence, winners, x.getKey()))
                .filter(x -> x.getValue() >= minEvidence)
                .collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

        if (unmatched.isEmpty())
            return Lists.newArrayList();

        List<Haplotype> haplotypes = unmatched.entrySet().stream()
                .map(x -> Haplotype.create(evidence.getAminoAcidIndexList(), new Pair(x.getKey(), x.getValue()), aminoAcidCount))
                .collect(Collectors.toList());

        return haplotypes;
    }
}
