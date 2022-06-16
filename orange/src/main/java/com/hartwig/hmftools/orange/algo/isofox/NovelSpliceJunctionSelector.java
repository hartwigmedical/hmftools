package com.hartwig.hmftools.orange.algo.isofox;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.orange.algo.linx.DNAFusionEvaluator;

import org.jetbrains.annotations.NotNull;

final class NovelSpliceJunctionSelector {

    private NovelSpliceJunctionSelector() {
    }

    @NotNull
    public static List<NovelSpliceJunction> selectSkippedExons(@NotNull List<NovelSpliceJunction> junctions,
            @NotNull List<LinxFusion> linxFusions, @NotNull KnownFusionCache knownFusionCache) {
        List<NovelSpliceJunction> result = Lists.newArrayList();

        for (NovelSpliceJunction junction : junctions) {
            if (knownFusionCache.hasExonDelDup(junction.geneName())) {
                boolean isTypeMatch = junction.type().equals("SKIPPED_EXONS");
                boolean hasSufficientFragments = junction.fragmentCount() > 5;
                boolean hasLimitedCohortFreq = junction.cohortFrequency() < 30;
                boolean hasReportedLinxFusion = DNAFusionEvaluator.hasFusion(linxFusions, junction.geneName(), junction.geneName());
                if (isTypeMatch && hasSufficientFragments && hasLimitedCohortFreq && !hasReportedLinxFusion) {
                    result.add(junction);
                }
            }
        }

        return result;
    }

    @NotNull
    public static List<NovelSpliceJunction> selectNovelExonsIntrons(@NotNull List<NovelSpliceJunction> junctions,
            @NotNull List<DriverGene> driverGenes) {
        List<NovelSpliceJunction> result = Lists.newArrayList();

        Set<String> drivers = Sets.newHashSet();
        for (DriverGene driverGene : driverGenes) {
            drivers.add(driverGene.gene());
        }

        for (NovelSpliceJunction junction : junctions) {
            if (drivers.contains(junction.geneName())) {
                boolean isTypeMatch = junction.type().equals("NOVEL_INTRON") || junction.type().equals("NOVEL_EXON");
                boolean hasSufficientFragments = junction.fragmentCount() > 5;
                boolean hasLimitedCohortFreq = junction.cohortFrequency() < 10;
                if (isTypeMatch && hasSufficientFragments && hasLimitedCohortFreq) {
                    result.add(junction);
                }
            }
        }

        return result;
    }
}
