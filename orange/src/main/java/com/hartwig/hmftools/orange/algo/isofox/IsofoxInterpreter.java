package com.hartwig.hmftools.orange.algo.isofox;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;

public final class IsofoxInterpreter {

    private IsofoxInterpreter() {
    }

    @NotNull
    public static IsofoxInterpretedData interpret(@NotNull IsofoxData isofox, @NotNull LinxData linx, @NotNull List<DriverGene> driverGenes,
            @NotNull KnownFusionCache knownFusionCache) {
        return ImmutableIsofoxInterpretedData.builder()
                .summary(isofox.summary())
                .allGeneExpressions(isofox.geneExpressions())
                .reportableHighExpression(selectHighExpressionGenes(isofox.geneExpressions(), driverGenes))
                .reportableLowExpression(selectLowExpressionGenes(isofox.geneExpressions(), driverGenes))
                .allFusions(isofox.fusions())
                .reportableNovelKnownFusions(selectNovelKnownFusions(isofox.fusions(), linx.reportableFusions(), knownFusionCache))
                .reportableNovelPromiscuousFusions(selectNovelPromiscuousFusions(isofox.fusions(),
                        linx.reportableFusions(),
                        knownFusionCache))
                .allNovelSpliceJunctions(isofox.novelSpliceJunctions())
                .reportableSkippedExons(selectSkippedExons(isofox.novelSpliceJunctions(), linx.reportableFusions(), knownFusionCache))
                .reportableNovelExonsIntrons(selectNovelExonsIntrons(isofox.novelSpliceJunctions(), driverGenes))
                .build();
    }

    @NotNull
    private static List<GeneExpression> selectHighExpressionGenes(@NotNull List<GeneExpression> expressions,
            @NotNull List<DriverGene> driverGenes) {
        Set<String> oncogenes = extractGenesOfType(driverGenes, DriverCategory.ONCO);

        List<GeneExpression> result = Lists.newArrayList();
        for (GeneExpression expression : expressions) {
            if (oncogenes.contains(expression.geneName()) && expression.percentileCohort() > 0.9 && expression.percentileCancer() > 0.9) {
                result.add(expression);
            }
        }

        return result;
    }

    @NotNull
    private static List<GeneExpression> selectLowExpressionGenes(@NotNull List<GeneExpression> expressions,
            @NotNull List<DriverGene> driverGenes) {
        Set<String> suppressors = extractGenesOfType(driverGenes, DriverCategory.TSG);

        List<GeneExpression> result = Lists.newArrayList();
        for (GeneExpression expression : expressions) {
            if (suppressors.contains(expression.geneName()) && expression.percentileCohort() < 0.05
                    && expression.percentileCancer() < 0.05) {
                result.add(expression);
            }
        }

        return result;
    }

    @NotNull
    private static List<RnaFusion> selectNovelKnownFusions(@NotNull List<RnaFusion> rnaFusions, @NotNull List<LinxFusion> linxFusions,
            @NotNull KnownFusionCache knownFusionCache) {
        List<RnaFusion> result = Lists.newArrayList();

        return result;
    }

    @NotNull
    private static List<RnaFusion> selectNovelPromiscuousFusions(@NotNull List<RnaFusion> rnaFusions, @NotNull List<LinxFusion> linxFusions,
            @NotNull KnownFusionCache knownFusionCache) {
        List<RnaFusion> result = Lists.newArrayList();

        return result;
    }

    @NotNull
    private static List<NovelSpliceJunction> selectSkippedExons(@NotNull List<NovelSpliceJunction> novelSpliceJunctions,
            @NotNull List<LinxFusion> linxFusions, @NotNull KnownFusionCache knownFusionCache) {
        List<NovelSpliceJunction> result = Lists.newArrayList();

        return result;
    }

    @NotNull
    private static List<NovelSpliceJunction> selectNovelExonsIntrons(@NotNull List<NovelSpliceJunction> novelSpliceJunctions,
            @NotNull List<DriverGene> driverGenes) {
        List<NovelSpliceJunction> result = Lists.newArrayList();

        return result;
    }

    @NotNull
    private static Set<String> extractGenesOfType(@NotNull List<DriverGene> driverGenes, @NotNull DriverCategory categoryToInclude) {
        Set<String> filtered = Sets.newHashSet();
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.likelihoodType() == categoryToInclude) {
                filtered.add(driverGene.gene());
            }
        }
        return filtered;
    }
}
