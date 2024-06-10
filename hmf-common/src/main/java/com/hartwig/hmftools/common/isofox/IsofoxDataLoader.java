package com.hartwig.hmftools.common.isofox;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaStatistics;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class IsofoxDataLoader
{
    private static final Logger LOGGER = LogManager.getLogger(IsofoxDataLoader.class);

    public static IsofoxData load(final String isofoxCancerType, final String isofoxGeneDistributionCsv,
            final String isofoxAltSjCohortCsv, final String isofoxSummaryCsv, final String isofoxGeneDataCsv,
            final String isofoxFusionCsv, final String isofoxAltSpliceJunctionCsv) throws IOException
    {
        LOGGER.info("Loading ISOFOX data from {} using cancer type '{}'", new File(isofoxSummaryCsv).getParent(), isofoxCancerType);

        List<String> summaryLines = Files.readAllLines(Paths.get(isofoxSummaryCsv));
        RnaStatistics summary = RnaStatistics.fromLines(summaryLines);

        LOGGER.info((" Loaded summary from " + isofoxSummaryCsv));

        GeneExpressionDistributionData geneDistributionData = new GeneExpressionDistributionData(isofoxGeneDistributionCsv);
        if(isofoxCancerType != null && !geneDistributionData.configuredCancerTypes().contains(isofoxCancerType))
        {
            throw new IllegalStateException("Cancer type does not exist as cohort in isofox gene distribution data: " + isofoxCancerType);
        }

        List<GeneExpression> geneExpressions =
                GeneExpressionLoader.loadGeneExpression(isofoxGeneDataCsv, geneDistributionData, isofoxCancerType);
        LOGGER.info(" Loaded {} gene expressions from {}", geneExpressions.size(), isofoxGeneDataCsv);

        List<RnaFusion> fusions = IsofoxFusionLoader.load(isofoxFusionCsv);
        LOGGER.info(" Loaded {} fusions from {}", fusions.size(), isofoxFusionCsv);

        AltSjCohortData altSjData = new AltSjCohortData(isofoxAltSjCohortCsv);
        List<NovelSpliceJunction> novelSpliceJunctions = NovelSpliceJunctionLoader.load(isofoxAltSpliceJunctionCsv, altSjData);
        LOGGER.info(" Loaded {} novel splice junctions from {}", novelSpliceJunctions.size(), isofoxAltSpliceJunctionCsv);

        return ImmutableIsofoxData.builder()
                .summary(summary)
                .geneExpressions(geneExpressions)
                .fusions(fusions)
                .novelSpliceJunctions(novelSpliceJunctions)
                .build();
    }
}
