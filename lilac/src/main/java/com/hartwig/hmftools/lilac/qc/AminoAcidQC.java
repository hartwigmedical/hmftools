package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.LOG_UNMATCHED_HAPLOTYPE_SUPPORT;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_UNMATCHED_HAPLOTYPE_SUPPORT;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.SequenceCount;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class AminoAcidQC
{
    public final int UnusedAminoAcids;
    public final int UnusedAminoAcidMaxSupport;

    public AminoAcidQC(int unusedAminoAcids, int unusedAminoAcidMaxSupport)
    {
        UnusedAminoAcids = unusedAminoAcids;
        UnusedAminoAcidMaxSupport = unusedAminoAcidMaxSupport;
    }

    public final List<String> header()
    {
        return Lists.newArrayList("unusedAminoAcids", "unusedAminoAcidMaxSupport");
    }

    public final List<String> body()
    {
        return Lists.newArrayList(String.valueOf(UnusedAminoAcids), String.valueOf(UnusedAminoAcidMaxSupport));
    }

    public static AminoAcidQC create(
            final List<HlaSequenceLoci> winners, final SequenceCount aminoAcidCount, final List<Haplotype> unmatchedHaplotypes)
    {
        int unused = 0;
        int largest = 0;

        for(Integer locus : aminoAcidCount.heterozygousLoci())
        {
            // ignore amino acids in any unmatched haplotype
            if(unmatchedHaplotypes.stream().anyMatch(x -> positionWithin(locus, x.StartLocus, x.EndLocus)))
                continue;

            Map<String,Integer> expected = aminoAcidCount.get(locus);

            Set<String> actualSequences =  winners.stream()
                    .filter(x -> locus < x.getSequences().size()).map(x -> x.sequence(locus)).collect(Collectors.toSet());

            if(actualSequences.stream().anyMatch(x -> x.equals(WILD_STR)))
                continue;

            for(Map.Entry<String,Integer> entry : expected.entrySet())
            {
                if(actualSequences.contains(entry.getKey()))
                    continue;

                int count = entry.getValue();

                if(count < LOG_UNMATCHED_HAPLOTYPE_SUPPORT)
                    continue;

                LL_LOGGER.warn("  UNMATCHED_AMINO_ACID - amino acid sequence({} count={} locus={}) not in winning solution",
                        entry.getKey(), count, locus);

                if(count >= WARN_UNMATCHED_HAPLOTYPE_SUPPORT)
                {
                    ++unused;
                    largest = max(largest, count);
                }
            }
        }

        return new AminoAcidQC(unused, largest);
    }
}
