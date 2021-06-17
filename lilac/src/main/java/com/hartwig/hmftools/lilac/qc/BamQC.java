package com.hartwig.hmftools.lilac.qc;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_LOW_COVERAGE_THRESHOLD;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.read.BamReader;
import com.hartwig.hmftools.lilac.read.BamRecordReader;
import com.hartwig.hmftools.lilac.read.Indel;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.jetbrains.annotations.NotNull;

public class BamQC
{
    public final int DiscardedAlignmentFragments;
    public final int DiscardedIndels;
    public final int DiscardedIndelMaxFrags;

    private final Map<String,int[]> mGeneBaseDepth;

    private static final int MIN_SUPPORT = 3;

    public BamQC(int discardedAlignmentFragments, int discardedIndels, int discardedIndelMaxFrags, final Map<String,int[]> geneBaseDepth)
    {
        DiscardedAlignmentFragments = discardedAlignmentFragments;
        DiscardedIndels = discardedIndels;
        DiscardedIndelMaxFrags = discardedIndelMaxFrags;

        mGeneBaseDepth = geneBaseDepth;
    }

    public List<String> header()
    {
        return Lists.newArrayList("DiscardedIndels", "DiscardedIndelMaxFrags", "DiscardedAlignmentFragments",
                "A_LowCoverageBases", "B_LowCoverageBases", "C_LowCoverageBases");
    }

    @NotNull
    public List<String> body()
    {
        int aLowCoverage = (int)Arrays.stream(mGeneBaseDepth.get(HLA_A)).filter(x -> x < WARN_LOW_COVERAGE_THRESHOLD).count();
        int bLowCoverage = (int)Arrays.stream(mGeneBaseDepth.get(HLA_B)).filter(x -> x < WARN_LOW_COVERAGE_THRESHOLD).count();
        int cLowCoverage = (int)Arrays.stream(mGeneBaseDepth.get(HLA_C)).filter(x -> x < WARN_LOW_COVERAGE_THRESHOLD).count();

        return Lists.newArrayList(
                String.valueOf(DiscardedIndels), String.valueOf(DiscardedIndelMaxFrags), String.valueOf(DiscardedAlignmentFragments),
                String.valueOf(aLowCoverage), String.valueOf(bLowCoverage), String.valueOf(cLowCoverage));
    }

    public static BamQC create(final BamReader reader, final Map<String,int[]> geneBaseDepth)
    {
        Map<Indel,Integer> fragmentsWithUnmatchedPonIndel = reader.unmatchedPonIndels(MIN_SUPPORT);
        Map<Indel,Integer> fragmentsWithUnmatchedIndel = reader.unmatchedIndels(MIN_SUPPORT);

        fragmentsWithUnmatchedIndel.entrySet().forEach(x -> LL_LOGGER.warn(
                "  UNMATCHED_INDEL - {} fragments excluded with unmatched indel {}", x.getValue(), x.getKey().toString()));

        fragmentsWithUnmatchedPonIndel.entrySet().forEach(x -> LL_LOGGER.warn(
                "  UNMATCHED_PON_INDEL - {} fragments excluded with unmatched PON indel {}", x.getValue(), x.getKey().toString()));

        return new BamQC(
                reader.alignmentFiltered(),
                fragmentsWithUnmatchedIndel.size(),
                fragmentsWithUnmatchedIndel.values().stream().mapToInt(x -> x).max().orElse(0),
                geneBaseDepth);
    }

}
