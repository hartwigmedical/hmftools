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
    private final int mDiscardedAlignmentFragments;
    private final int mDiscardedIndelFragments;
    private final int mDiscardedIndelMaxCount;

    private static final int MIN_SUPPORT = 3;

    public BamQC(
            int discardedAlignmentFragments, int discardedIndelFragments, int discardedIndelMaxCount)
    {
        mDiscardedAlignmentFragments = discardedAlignmentFragments;
        mDiscardedIndelFragments = discardedIndelFragments;
        mDiscardedIndelMaxCount = discardedIndelMaxCount;
    }

    public final List<String> header()
    {
        return Lists.newArrayList("discardedIndelFragments", "discardedIndelMaxCount", "discardedAlignmentFragments");
    }

    @NotNull
    public final List<String> body()
    {
        return Lists.newArrayList(
                String.valueOf(mDiscardedIndelFragments),
                String.valueOf(mDiscardedIndelMaxCount),
                String.valueOf(mDiscardedAlignmentFragments));
    }

    public final int getDiscardedIndelFragments()
    {
        return mDiscardedIndelFragments;
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
