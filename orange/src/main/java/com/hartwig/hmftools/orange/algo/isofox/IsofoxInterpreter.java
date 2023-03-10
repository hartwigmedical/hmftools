package com.hartwig.hmftools.orange.algo.isofox;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.util.List;

public class IsofoxInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(IsofoxInterpreter.class);

    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;
    @NotNull
    private final LinxRecord linx;

    public IsofoxInterpreter(@NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache,
                             @NotNull final LinxRecord linx) {
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
        this.linx = linx;
    }

    @NotNull
    public IsofoxInterpretedData interpret(@NotNull IsofoxData isofox) {
        List<GeneExpression> highExpressionGenes = ExpressionSelector.selectHighExpressionGenes(isofox.geneExpressions(), driverGenes);
        LOGGER.info(" Found {} genes with high expression", highExpressionGenes.size());

        List<GeneExpression> lowExpressionGenes = ExpressionSelector.selectLowExpressionGenes(isofox.geneExpressions(), driverGenes);
        LOGGER.info(" Found {} genes with low expression", lowExpressionGenes.size());

        List<RnaFusion> novelKnownFusions =
                RNAFusionSelector.selectNovelKnownFusions(isofox.fusions(), linx.allSomaticFusions(), knownFusionCache);
        LOGGER.info(" Found {} novel known fusions in RNA", novelKnownFusions.size());

        List<RnaFusion> novelPromiscuousFusions =
                RNAFusionSelector.selectNovelPromiscuousFusions(isofox.fusions(), linx.allSomaticFusions(), knownFusionCache);
        LOGGER.info(" Found {} novel promiscuous fusions in RNA", novelPromiscuousFusions.size());

        List<NovelSpliceJunction> suspiciousSkippedExons =
                NovelSpliceJunctionSelector.selectSkippedExons(isofox.novelSpliceJunctions(), linx.allSomaticFusions(), knownFusionCache);
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
