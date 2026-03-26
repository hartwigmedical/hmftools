package com.hartwig.hmftools.orange.algo.isofox;

import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.geneNames;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.ALT_SJ_MIN_FRAGMENTS;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.HIGH_EXPRESSION_PERCENTILE_CUTOFF;
import static com.hartwig.hmftools.orange.algo.OrangeConstants.LOW_EXPRESSION_PERCENTILE_CUTOFF;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.orange.conversion.IsofoxConversion;

public class IsofoxInterpreter
{
    private final LinxRecord mLinxRecord;

    public IsofoxInterpreter(final LinxRecord linx)
    {
        mLinxRecord = linx;
    }

    public IsofoxRecord interpret(final IsofoxData isofox)
    {
        List<GeneExpression> highExpressionGenes = Lists.newArrayList();
        List<GeneExpression> lowExpressionGenes = Lists.newArrayList();

        findExpressionOutliers(isofox.geneExpressions(), highExpressionGenes, lowExpressionGenes);

        List<RnaFusion> fusions = findFusions(isofox.fusions(), mLinxRecord.fusions());

        List<NovelSpliceJunction> novelSpliceJunctions = findNovelSplceJunctions(isofox.novelSpliceJunctions(), fusions);

        return ImmutableIsofoxRecord.builder()
                .summary(IsofoxConversion.convert(isofox.summary()))
                .highExpressionGenes(highExpressionGenes)
                .lowExpressionGenes(lowExpressionGenes)
                .fusions(fusions)
                .novelSpliceJunctions(novelSpliceJunctions)
                .build();
    }

    @VisibleForTesting
    public static void findExpressionOutliers(
            final List<com.hartwig.hmftools.common.rna.GeneExpression> geneExpressions,
            final List<GeneExpression> highExpressionGenes, final List<GeneExpression> lowExpressionGenes)
    {
        for(com.hartwig.hmftools.common.rna.GeneExpression geneExpression : geneExpressions)
        {
            if(geneExpression.medianTpmCohort() <= 0) // points to an issue with annotation in Isofox
                continue;

            if(geneExpression.reportedStatus() == ReportedStatus.REPORTED)
                highExpressionGenes.add(IsofoxConversion.convert(geneExpression));
        }
    }

    public static List<NovelSpliceJunction> findNovelSplceJunctions(
            final List<com.hartwig.hmftools.common.rna.NovelSpliceJunction> altSpliceJunctions,
            final List<RnaFusion> fusions)
    {
        List<NovelSpliceJunction> spliceJunctions = Lists.newArrayList();

        for(com.hartwig.hmftools.common.rna.NovelSpliceJunction altSpliceJunction : altSpliceJunctions)
        {
            // all in the passing file should be shown
            spliceJunctions.add(IsofoxConversion.convert(altSpliceJunction));
        }

        return spliceJunctions;
    }

    @VisibleForTesting
    public static List<RnaFusion> findFusions(final List<com.hartwig.hmftools.common.rna.RnaFusion> rnaFusions, final List<LinxFusion> linxFusions)
    {
        List<RnaFusion> reportableFusions = Lists.newArrayList();

        for(com.hartwig.hmftools.common.rna.RnaFusion rnaFusion : rnaFusions)
        {
            String[] geneNames = geneNames(rnaFusion);

            if(geneNames == null || geneNames.length != 2 || geneNames[0].isEmpty()|| geneNames[1].isEmpty())
                continue;

            String geneUp = geneNames[0];
            String geneDown = geneNames[1];

            if(matchesDnaFusion(linxFusions, geneUp, geneDown))
                continue;

            // can expand the set of interesting fusions once Isofox sets this field similarly to Linx and a likelihood
            if(rnaFusion.knownType() == KNOWN_PAIR)
            {
                reportableFusions.add(IsofoxConversion.convert(rnaFusion));
            }

            /*
            if(hasKnownPairGene(rnaFusion.knownType()))
            {
                reportableFusions.add(IsofoxConversion.convert(rnaFusion));
            }
            else if(hasPromiscousGene(rnaFusion.knownType()))
            {
                boolean validForPromiscuous = false;

                if(rnaFusion.svType() == BND || rnaFusion.svType() == INV)
                    validForPromiscuous = true;
                else if(rnaFusion.svType() == DEL || rnaFusion.svType() == DUP)
                    validForPromiscuous = abs(rnaFusion.positionUp() - rnaFusion.positionDown()) > RNA_MIN_PROMISCOUS_DISTANCE;

                if(validForPromiscuous)
                    reportableFusions.add(IsofoxConversion.convert(rnaFusion));
            }
            */
        }

        return reportableFusions;
    }

    private static boolean matchesDnaFusion(final List<LinxFusion> linxFusions, final String geneUp, final String geneDown)
    {
        for(LinxFusion linxFusion : linxFusions)
        {
            if(linxFusion.geneUp().equals(geneUp) && linxFusion.geneDown().equals(geneDown))
            {
                return true;
            }
        }
        return false;
    }
}
