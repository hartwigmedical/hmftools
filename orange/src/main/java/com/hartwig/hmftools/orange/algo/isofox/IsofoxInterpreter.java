package com.hartwig.hmftools.orange.algo.isofox;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class IsofoxInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(IsofoxInterpreter.class);

    private IsofoxInterpreter() {
    }

    @NotNull
    public static IsofoxInterpretedData interpret(@NotNull IsofoxData isofox, @NotNull LinxInterpretedData linx,
            @NotNull List<DriverGene> driverGenes, @NotNull KnownFusionCache knownFusionCache) {
        List<LinxFusion> allFusions = Lists.newArrayList();
        allFusions.addAll(linx.reportableFusions());
        allFusions.addAll(linx.allFusions());

        List<GeneExpression> highExpressionGenes = ExpressionSelector.selectHighExpressionGenes(isofox.geneExpressions(), driverGenes);
        LOGGER.info(" Found {} genes with high expression", highExpressionGenes.size());

        List<GeneExpression> lowExpressionGenes = ExpressionSelector.selectLowExpressionGenes(isofox.geneExpressions(), driverGenes);
        LOGGER.info(" Found {} genes with low expression", lowExpressionGenes.size());

        List<RnaFusion> novelKnownFusions = RNAFusionSelector.selectNovelKnownFusions(isofox.fusions(), allFusions, knownFusionCache);
        LOGGER.info(" Found {} novel known fusions in RNA", novelKnownFusions.size());

        List<RnaFusion> novelPromiscuousFusions =
                RNAFusionSelector.selectNovelPromiscuousFusions(isofox.fusions(), allFusions, knownFusionCache);
        LOGGER.info(" Found {} novel promiscuous fusions in RNA", novelPromiscuousFusions.size());

        List<NovelSpliceJunction> suspiciousSkippedExons =
                NovelSpliceJunctionSelector.selectSkippedExons(isofox.novelSpliceJunctions(), allFusions, knownFusionCache);
        LOGGER.info(" Found {} suspicious skipped exons in RNA", suspiciousSkippedExons.size());

        List<NovelSpliceJunction> suspiciousNovelExonsIntrons =
                NovelSpliceJunctionSelector.selectNovelExonsIntrons(isofox.novelSpliceJunctions(), driverGenes);
        LOGGER.info(" Found {} suspicious novel exons/introns in RNA", suspiciousNovelExonsIntrons.size());

        return ImmutableIsofoxInterpretedData.builder()
                .summary(isofox.summary())
                .allGeneExpressions(isofox.geneExpressions())
                .reportableHighExpression(highExpressionGenes)
                .reportableLowExpression(lowExpressionGenes)
                .allFusions(isofox.fusions())
                .reportableNovelKnownFusions(novelKnownFusions)
                .reportableNovelPromiscuousFusions(novelPromiscuousFusions)
                .allNovelSpliceJunctions(isofox.novelSpliceJunctions())
                .reportableSkippedExons(suspiciousSkippedExons)
                .reportableNovelExonsIntrons(suspiciousNovelExonsIntrons)
                .build();
    }
}
