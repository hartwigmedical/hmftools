package com.hartwig.hmftools.orange.algo.isofox;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.common.rna.RnaStatisticFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class IsofoxDataLoader
{
    private static final Logger LOGGER = LogManager.getLogger(IsofoxDataLoader.class);

    public static IsofoxData load(
            final String isofoxCancerType, final String isofoxSummaryCsv, final String isofoxGeneDataCsv,
            final String isofoxFusionCsv, final String isofoxAltSpliceJunctionCsv) throws IOException
    {
        LOGGER.info("Loading ISOFOX data from {}", new File(isofoxSummaryCsv).getParent(), isofoxCancerType);

        List<String> summaryLines = Files.readAllLines(Paths.get(isofoxSummaryCsv));
        RnaStatistics summary = RnaStatisticFile.fromLines(summaryLines);

        LOGGER.info((" Loaded summary from " + isofoxSummaryCsv));

        List<GeneExpression> geneExpressions = GeneExpressionFile.read(isofoxGeneDataCsv);

        LOGGER.info(" Loaded {} gene expressions from {}", geneExpressions.size(), isofoxGeneDataCsv);

        List<RnaFusion> fusions = RnaFusionFile.read(isofoxFusionCsv);
        LOGGER.info(" Loaded {} fusions from {}", fusions.size(), isofoxFusionCsv);

        List<NovelSpliceJunction> novelSpliceJunctions = NovelSpliceJunctionFile.read(isofoxAltSpliceJunctionCsv);
        LOGGER.info(" Loaded {} novel splice junctions from {}", novelSpliceJunctions.size(), isofoxAltSpliceJunctionCsv);

        return ImmutableIsofoxData.builder()
                .summary(summary)
                .geneExpressions(geneExpressions)
                .fusions(fusions)
                .novelSpliceJunctions(novelSpliceJunctions)
                .build();
    }
}
