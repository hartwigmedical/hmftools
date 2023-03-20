package com.hartwig.hmftools.orange.algo.isofox;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxInterpretedData;
import com.hartwig.hmftools.datamodel.isofox.IsofoxInterpretedData;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.rna.GeneExpression;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.IsofoxConversion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IsofoxInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(IsofoxInterpreter.class);

    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;
    @NotNull
    private final LinxRecord linx;

    public IsofoxInterpreter(@NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache,
            @NotNull final LinxRecord linx)
    {
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
        this.linx = linx;
    }

    @NotNull
    public IsofoxInterpretedData interpret(@NotNull IsofoxData isofox) {
        List<GeneExpression> geneExpressions = ConversionUtil.mapToList(isofox.geneExpressions(), IsofoxConversion::convert);
        List<GeneExpression> highExpressionGenes = ExpressionSelector.selectHighExpressionGenes(geneExpressions, driverGenes);
        LOGGER.info(" Found {} genes with high expression", highExpressionGenes.size());

        List<GeneExpression> lowExpressionGenes = ExpressionSelector.selectLowExpressionGenes(geneExpressions, driverGenes);
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
                .summary(IsofoxConversion.convert(isofox.summary()))
                .allGeneExpressions(geneExpressions)
                .reportableHighExpression(highExpressionGenes)
                .reportableLowExpression(lowExpressionGenes)
                .allFusions(ConversionUtil.mapToIterable(isofox.fusions(), IsofoxConversion::convert))
                .reportableNovelKnownFusions(ConversionUtil.mapToIterable(novelKnownFusions, IsofoxConversion::convert))
                .reportableNovelPromiscuousFusions(ConversionUtil.mapToIterable(novelPromiscuousFusions, IsofoxConversion::convert))
                .allNovelSpliceJunctions(ConversionUtil.mapToIterable(isofox.novelSpliceJunctions(), IsofoxConversion::convert))
                .reportableSkippedExons(ConversionUtil.mapToIterable(suspiciousSkippedExons, IsofoxConversion::convert))
                .reportableNovelExonsIntrons(ConversionUtil.mapToIterable(suspiciousNovelExonsIntrons, IsofoxConversion::convert))
                .build();
    }
}
