package com.hartwig.hmftools.lilac.qc;

import static java.lang.String.format;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_LOW_COVERAGE_DEPTH;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.BamReader;
import com.hartwig.hmftools.lilac.read.Indel;

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

    public int totalLowCoverage() { return GeneLowCoverageCounts.values().stream().mapToInt(Integer::intValue).sum(); }

    public List<String> header()
    {
        return List.of("DiscardedIndels", "DiscardedIndelMaxFrags", "DiscardedAlignmentFragments", "LowCoverageBases");
    }

    public List<String> body()
    {
        StringJoiner lowCoverageBasesStrBuilder = new StringJoiner(";");
        GENE_CACHE.GeneNames.stream()
                .filter(x -> !x.isPseudo())
                .map(gene -> format("%s=%d", gene.shortName(), GeneLowCoverageCounts.getOrDefault(gene, 0)))
                .forEach(lowCoverageBasesStrBuilder::add);

        return List.of(String.valueOf(DiscardedIndels), String.valueOf(DiscardedIndelMaxFrags), String.valueOf(DiscardedAlignmentFragments), lowCoverageBasesStrBuilder.toString());
    }

    public static BamQC create(final BamReader reader, final Map<HlaGene, int[]> geneBaseDepth)
    {
        Map<Indel, Integer> fragmentsWithUnmatchedPonIndel = reader.unmatchedPonIndels(MIN_SUPPORT);
        Map<Indel, Integer> fragmentsWithUnmatchedIndel = reader.unmatchedIndels(MIN_SUPPORT);

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
