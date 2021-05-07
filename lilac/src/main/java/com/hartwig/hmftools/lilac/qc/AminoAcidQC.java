package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.max;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilackt.SequenceCount;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class AminoAcidQC
{
    private final int mUnusedAminoAcids;
    private final int mUnusedAminoAcidMaxSupport;
    private static final int MIN_COUNT = 3;

    public AminoAcidQC(int unusedAminoAcids, int unusedAminoAcidMaxSupport)
    {
        mUnusedAminoAcids = unusedAminoAcids;
        mUnusedAminoAcidMaxSupport = unusedAminoAcidMaxSupport;
    }

    public final List<String> header()
    {
        return Lists.newArrayList("unusedAminoAcids", "unusedAminoAcidMaxSupport");
    }

    public final List<String> body()
    {
        return Lists.newArrayList(String.valueOf(mUnusedAminoAcids), String.valueOf(mUnusedAminoAcidMaxSupport));
    }

    public static AminoAcidQC create(final Set<HlaSequenceLoci> winners, final SequenceCount aminoAcidCount)
    {
        int unused = 0;
        int largest = 0;

        for (Integer locus : aminoAcidCount.heterozygousLoci())
        {
            Map<String,Integer> expected = aminoAcidCount.get(locus);
            Set<String> actualSequences =  winners.stream()
                    .filter(x -> locus < x.getSequences().size()).map(x -> x.sequence(locus)).collect(Collectors.toSet());

            // CHECK
            // val actual = winners.filter { locus < it.sequences.size }.map { it.sequence(locus) }.toSet()
            for(Map.Entry<String,Integer> entry : expected.entrySet())
            {
                int count = entry.getValue();
                if(count >= MIN_COUNT && !actualSequences.contains(entry.getKey()))
                {
                    ++unused;
                    largest = max(largest, count);
                    LL_LOGGER.warn("  UNMATCHED_AMINO_ACID - amino acid sequence({})' with count({}) support at locus({}) not in winning solution",
                            entry.getKey(), count, locus);
                }
            }
        }

        return new AminoAcidQC(unused, largest);
    }
}
