package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.LOG_UNMATCHED_HAPLOTYPE_SUPPORT;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceCount;

public class AminoAcidQC
{
    public final int UnusedAminoAcids;
    public final int UnusedAminoAcidMaxFrags;

    public AminoAcidQC(int unusedAminoAcids, int unusedAminoAcidMaxFrags)
    {
        UnusedAminoAcids = unusedAminoAcids;
        UnusedAminoAcidMaxFrags = unusedAminoAcidMaxFrags;
    }

    public static List<String> header()
    {
        return Lists.newArrayList("UnusedAminoAcids", "UnusedAminoAcidMaxFrags");
    }

    public List<String> body()
    {
        return Lists.newArrayList(String.valueOf(UnusedAminoAcids), String.valueOf(UnusedAminoAcidMaxFrags));
    }

    public static AminoAcidQC create(
            final Collection<HlaSequenceLoci> winners, final Collection<HlaSequenceLoci> hlaYSequenceLoci,
            final SequenceCount aminoAcidCount, final Collection<Haplotype> unmatchedHaplotypes)
    {
        int unused = 0;
        int largest = 0;

        for(Integer locus : aminoAcidCount.heterozygousLoci())
        {
            // ignore amino acids in any unmatched haplotype
            if(unmatchedHaplotypes.stream().anyMatch(x -> positionWithin(locus, x.StartLocus, x.EndLocus)))
                continue;

            // or those matching the winning or HLA-Y sequences

            Multiset<String> expected = aminoAcidCount.get(locus);

            Set<String> actualSequences = winners.stream()
                    .filter(x -> locus < x.getSequences().size()).map(x -> x.sequence(locus)).collect(Collectors.toSet());

            Set<String> hlaYSequences = hlaYSequenceLoci == null ? null : hlaYSequenceLoci.stream()
                    .filter(x -> locus < x.getSequences().size()).map(x -> x.sequence(locus)).collect(Collectors.toSet());

            if(actualSequences.stream().anyMatch(x -> x.equals(WILD_STR)))
                continue;

            for(Multiset.Entry<String> entry : expected.entrySet())
            {
                String aminoAcid = entry.getElement();

                if(actualSequences.contains(aminoAcid))
                    continue;

                if(hlaYSequences != null && hlaYSequences.contains(aminoAcid))
                    continue;

                int count = entry.getCount();

                if(count < LOG_UNMATCHED_HAPLOTYPE_SUPPORT)
                    continue;

                LL_LOGGER.warn("  UNMATCHED_AMINO_ACID - amino acid sequence({} count={} locus={}) not in winning solution",
                        aminoAcid, count, locus);

                ++unused;
                largest = max(largest, count);
            }
        }

        return new AminoAcidQC(unused, largest);
    }
}
