package com.hartwig.hmftools.lilac.qc;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.read.SAMRecordReader;
import com.hartwig.hmftools.lilac.read.Indel;

import java.util.List;
import java.util.Map;

import org.jetbrains.annotations.NotNull;

public class BamQC
{
    public final int DiscardedAlignmentFragments;
    public final int DiscardedIndelFragments;
    public final int DiscardedIndelMaxCount;

    private static final int MIN_SUPPORT = 3;

    public BamQC(
            int discardedAlignmentFragments, int discardedIndelFragments, int discardedIndelMaxCount)
    {
        DiscardedAlignmentFragments = discardedAlignmentFragments;
        DiscardedIndelFragments = discardedIndelFragments;
        DiscardedIndelMaxCount = discardedIndelMaxCount;
    }

    public final List<String> header()
    {
        return Lists.newArrayList("DiscardedIndelFragments", "DiscardedIndelMaxCount", "DiscardedAlignmentFragments");
    }

    @NotNull
    public final List<String> body()
    {
        return Lists.newArrayList(
                String.valueOf(DiscardedIndelFragments),
                String.valueOf(DiscardedIndelMaxCount),
                String.valueOf(DiscardedAlignmentFragments));
    }

    public final int getDiscardedIndelFragments()
    {
        return DiscardedIndelFragments;
    }

    public static BamQC create(final SAMRecordReader reader)
    {
        Map<Indel,Integer> fragmentsWithUnmatchedPonIndel = reader.unmatchedPonIndels(MIN_SUPPORT);
        Map<Indel,Integer> fragmentsWithUnmatchedIndel = reader.unmatchedIndels(MIN_SUPPORT);

        fragmentsWithUnmatchedIndel.entrySet().forEach(x -> LL_LOGGER.warn(
                "  UNMATCHED_INDEL - {} fragments excluded with unmatched indel {}", x.getValue(), x.getKey().toString()));

        fragmentsWithUnmatchedPonIndel.entrySet().forEach(x -> LL_LOGGER.warn(
                "  UNMATCHED_PON_INDEL - {} fragments excluded with unmatched PON indel {}", x.getValue(), x.getKey().toString()));

        return new BamQC(reader.alignmentFiltered(),
                fragmentsWithUnmatchedIndel.size(), fragmentsWithUnmatchedIndel.values().stream().mapToInt(x -> x).max().orElse(0));
    }
}
