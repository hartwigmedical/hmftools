package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.IsofoxConversion;

public class IsofoxInterpreter
{
    private final List<DriverGene> mDriverGenes;
    private final LinxRecord mLinxRecord;

    public IsofoxInterpreter(final List<DriverGene> driverGenes, final LinxRecord linx)
    {
        this.mDriverGenes = driverGenes;
        this.mLinxRecord = linx;
    }

    public IsofoxRecord interpret(final IsofoxData isofox)
    {
        List<GeneExpression> geneExpressions = ConversionUtil.mapToList(isofox.geneExpressions(), IsofoxConversion::convert);
        List<GeneExpression> highExpressionGenes = ExpressionSelector.selectHighExpressionGenes(geneExpressions, mDriverGenes);
        LOGGER.info(" Found {} genes with high expression", highExpressionGenes.size());

        List<GeneExpression> lowExpressionGenes = ExpressionSelector.selectLowExpressionGenes(geneExpressions, mDriverGenes);
        LOGGER.info(" Found {} genes with low expression", lowExpressionGenes.size());

        List<RnaFusion> novelKnownFusions = RnaFusionSelector.selectNovelKnownFusions(isofox.fusions(), mLinxRecord.fusions());
        LOGGER.info(" Found {} novel known fusions in RNA", novelKnownFusions.size());

        List<RnaFusion> novelPromiscuousFusions = RnaFusionSelector.selectNovelPromiscuousFusions(isofox.fusions(), mLinxRecord.fusions());
        LOGGER.info(" Found {} novel promiscuous fusions in RNA", novelPromiscuousFusions.size());

        List<NovelSpliceJunction> suspiciousSkippedExons = NovelSpliceJunctionSelector.selectSkippedExons(
                isofox.fusions(), isofox.novelSpliceJunctions(), mLinxRecord.fusions());
        LOGGER.info(" Found {} suspicious skipped exons in RNA", suspiciousSkippedExons.size());

        List<NovelSpliceJunction> suspiciousNovelExonsIntrons =
                NovelSpliceJunctionSelector.selectNovelExonsIntrons(isofox.novelSpliceJunctions(), mDriverGenes);
        LOGGER.info(" Found {} suspicious novel exons/introns in RNA", suspiciousNovelExonsIntrons.size());

        return ImmutableIsofoxRecord.builder()
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
