package com.hartwig.hmftools.lilac.qc;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_LOW_COVERAGE_DEPTH;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_C;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.BamReader;
import com.hartwig.hmftools.lilac.read.Indel;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jetbrains.annotations.NotNull;

public class BamQC
{
    public final int DiscardedAlignmentFragments;
    public final int DiscardedIndels;
    public final int DiscardedIndelMaxFrags;

    public final Map<HlaGene, Integer> GeneLowCoverageCounts;

    private static final int MIN_SUPPORT = 3;

    public BamQC(int discardedAlignmentFragments, int discardedIndels, int discardedIndelMaxFrags, final Map<HlaGene, int[]> geneBaseDepth)
    {
        DiscardedAlignmentFragments = discardedAlignmentFragments;
        DiscardedIndels = discardedIndels;
        DiscardedIndelMaxFrags = discardedIndelMaxFrags;

        GeneLowCoverageCounts = Maps.newHashMap();

        for(HlaGene gene : GENE_CACHE.GeneNames)
        {
            int lowCoverageCount = (int) Arrays.stream(geneBaseDepth.get(gene)).filter(x -> x < WARN_LOW_COVERAGE_DEPTH).count();
            GeneLowCoverageCounts.put(gene, lowCoverageCount);
        }
    }

    public int totalLowCoverage() { return GeneLowCoverageCounts.values().stream().mapToInt(x -> x.intValue()).sum(); }

    public List<String> header()
    {
        return Lists.newArrayList("DiscardedIndels", "DiscardedIndelMaxFrags", "DiscardedAlignmentFragments",
                "A_LowCoverageBases", "B_LowCoverageBases", "C_LowCoverageBases");
    }

    public List<String> body()
    {
        return Lists.newArrayList(
                String.valueOf(DiscardedIndels), String.valueOf(DiscardedIndelMaxFrags), String.valueOf(DiscardedAlignmentFragments),
                String.valueOf(GeneLowCoverageCounts.get(HLA_A)), String.valueOf(GeneLowCoverageCounts.get(HLA_B)),
                String.valueOf(GeneLowCoverageCounts.get(HLA_C)));
    }

    public static BamQC create(final BamReader reader, final Map<HlaGene, int[]> geneBaseDepth)
    {
        Map<Indel,Integer> fragmentsWithUnmatchedPonIndel = reader.unmatchedPonIndels(MIN_SUPPORT);
        Map<Indel,Integer> fragmentsWithUnmatchedIndel = reader.unmatchedIndels(MIN_SUPPORT);

        fragmentsWithUnmatchedIndel.entrySet().forEach(x -> LL_LOGGER.warn(
                "  UNMATCHED_INDEL - {} fragments excluded with unmatched indel {}", x.getValue(), x.getKey().toString()));

        fragmentsWithUnmatchedPonIndel.entrySet().forEach(x -> LL_LOGGER.debug(
                "  UNMATCHED_PON_INDEL - {} fragments excluded with unmatched PON indel {}", x.getValue(), x.getKey().toString()));

        return new BamQC(
                reader.filteredReadCount(),
                fragmentsWithUnmatchedIndel.size(),
                fragmentsWithUnmatchedIndel.values().stream().mapToInt(x -> x).max().orElse(0),
                geneBaseDepth);
    }

}
