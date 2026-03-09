package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.orange.algo.OrangeConstants.HIGH_EXPRESSION_PERCENTILE_CUTOFF;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.LOW_EXPRESSION_PERCENTILE_CUTOFF;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.orange.conversion.IsofoxConversion;

public class IsofoxInterpreter
{
    private final Map<String,DriverGene> mDriverGenes;
    private final LinxRecord mLinxRecord;

    public IsofoxInterpreter(final Map<String,DriverGene> driverGenes, final LinxRecord linx)
    {
        mDriverGenes = driverGenes;
        mLinxRecord = linx;
    }

    public IsofoxRecord interpret(final IsofoxData isofox)
    {
        List<GeneExpression> highExpressionGenes = Lists.newArrayList();
        List<GeneExpression> lowExpressionGenes = Lists.newArrayList();

        findExpressionOutliers(isofox.geneExpressions(), mDriverGenes, highExpressionGenes, lowExpressionGenes);

        List<RnaFusion> fusion = findFusions(isofox.fusions());

        /*
        List<RnaFusion> novelKnownFusions = RnaFusionSelector.selectKnownFusions(isofox.fusions(), mLinxRecord.fusions());
        LOGGER.info(" Found {} novel known fusions in RNA", novelKnownFusions.size());

        List<RnaFusion> novelPromiscuousFusions = RnaFusionSelector.selectNovelPromiscuousFusions(isofox.fusions(), mLinxRecord.fusions());
        LOGGER.info(" Found {} novel promiscuous fusions in RNA", novelPromiscuousFusions.size());
        */

        List<NovelSpliceJunction> novelSpliceJunctions = findNovelSplceJunctions(isofox.novelSpliceJunctions(), mDriverGenes);

        return ImmutableIsofoxRecord.builder()
                .summary(IsofoxConversion.convert(isofox.summary()))
                .highExpressionGenes(highExpressionGenes)
                .lowExpressionGenes(lowExpressionGenes)
                .fusions(fusion)
                .novelSpliceJunctions(novelSpliceJunctions)
                .build();
    }

    @VisibleForTesting
    public static void findExpressionOutliers(
            final List<com.hartwig.hmftools.common.rna.GeneExpression> geneExpressions, final Map<String,DriverGene> driverGenes,
            final List<GeneExpression> highExpressionGenes, final List<GeneExpression> lowExpressionGenes)
    {
        for(com.hartwig.hmftools.common.rna.GeneExpression geneExpression : geneExpressions)
        {
            DriverGene driverGene = driverGenes.get(geneExpression.geneName());

            if(driverGene == null)
                continue;

            if(driverGene.reportHighExpression() && geneExpression.percentileCohort() > HIGH_EXPRESSION_PERCENTILE_CUTOFF)
                highExpressionGenes.add(IsofoxConversion.convert(geneExpression));
            else if(driverGene.reportLowExpression() && geneExpression.percentileCohort() < LOW_EXPRESSION_PERCENTILE_CUTOFF)
                lowExpressionGenes.add(IsofoxConversion.convert(geneExpression));
        }
    }

    public static List<NovelSpliceJunction> findNovelSplceJunctions(
            final List<com.hartwig.hmftools.common.rna.NovelSpliceJunction> altSpliceJunctions,
            final Map<String,DriverGene> driverGenes)
    {
            List<NovelSpliceJunction> spliceJunctions = Lists.newArrayList();

            for(com.hartwig.hmftools.common.rna.NovelSpliceJunction altSpliceJunction : altSpliceJunctions)
            {
                DriverGene driverGene = driverGenes.get(altSpliceJunction.geneName());

                if(driverGene == null || !driverGene.reportNovelSpliceJunctions())
                    continue;

                // TODO: consider other criteria

                spliceJunctions.add(IsofoxConversion.convert(altSpliceJunction));

            /*
            for(RnaFusion rnaFusion : rnaFusions)
            {
                String geneUp = geneUp(rnaFusion);
                String geneDown = geneDown(rnaFusion);
                if(geneUp != null && geneDown != null && geneUp.equals(geneDown)
                        && rnaFusion.knownType() == KnownFusionType.KNOWN_PAIR)
                {
                    for(NovelSpliceJunction junction : junctions)
                    {
                        if(!junction.geneName().equals(geneUp))
                            continue;
                        if(junction.type() != SKIPPED_EXONS)
                            continue;

                        boolean isTypeMatch = junction.type() == SKIPPED_EXONS;
                        boolean hasSufficientFragments = junction.fragmentCount() > 5;
                        boolean hasLimitedCohortFreq = junction.cohortFrequency() < 30;
                        boolean hasReportedLinxFusion = hasFusion(linxFusions, junction.geneName(), junction.geneName());
                        if(isTypeMatch && hasSufficientFragments && hasLimitedCohortFreq && !hasReportedLinxFusion)
                        {
                            result.add(junction);
                        }
                    }
                }
            }

            if(drivers.contains(junction.geneName()))
            {
                boolean isTypeMatch =
                        junction.type() == AltSpliceJunctionType.NOVEL_INTRON || junction.type() == AltSpliceJunctionType.NOVEL_EXON;
                boolean hasSufficientFragments = junction.fragmentCount() > 5;
                boolean hasLimitedCohortFreq = junction.cohortFrequency() < 10;
                if(isTypeMatch && hasSufficientFragments && hasLimitedCohortFreq)
                {
                    result.add(junction);
                }
            }

            */

        }

        return spliceJunctions;
    }

    @VisibleForTesting
    public static List<RnaFusion> findFusions(final List<com.hartwig.hmftools.common.rna.RnaFusion> fusions)
    {
        List<RnaFusion> reportableFusions = Lists.newArrayList();

        return reportableFusions;
    }
}
