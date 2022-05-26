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
import org.jetbrains.annotations.NotNull;

public final class IsofoxDataLoader {

    private static final Logger LOGGER = LogManager.getLogger(IsofoxDataLoader.class);

    private IsofoxDataLoader() {
    }

    @NotNull
    public static IsofoxData load(@NotNull String isofoxCancerType, @NotNull String isofoxGeneDistributionCsv,
            @NotNull String isofoxAltSjCohortCsv, @NotNull String isofoxSummaryCsv, @NotNull String isofoxGeneDataCsv,
            @NotNull String isofoxFusionCsv, @NotNull String isofoxAltSpliceJunctionCsv) throws IOException {
        LOGGER.info("Loading ISOFOX data from {} using cancer type '{}'", new File(isofoxSummaryCsv).getParent(), isofoxCancerType);

        List<String> summaryLines = Files.readAllLines(Paths.get(isofoxSummaryCsv));
        RnaStatistics summary = RnaStatistics.fromCsv(summaryLines.get(1));

        LOGGER.info((" Loaded summary from " + isofoxSummaryCsv));

        GeneExpressionDistributionData geneDistributionData = new GeneExpressionDistributionData(isofoxGeneDistributionCsv);
        List<GeneExpression> geneExpressions =
                GeneExpressionLoader.loadGeneExpression(isofoxGeneDataCsv, geneDistributionData, isofoxCancerType);
        LOGGER.info(" Loaded {} gene expressions from {}", geneExpressions.size(), isofoxGeneDataCsv);

        List<RnaFusion> fusions = IsofoxFusionLoader.load(isofoxFusionCsv);
        LOGGER.info(" Loaded {} fusions from {}", fusions.size(), isofoxFusionCsv);

        AltSjCohortData altSjData = new AltSjCohortData(isofoxAltSjCohortCsv);
        List<NovelSpliceJunction> novelSpliceJunctions = IsofoxNovelSpliceJunctionLoader.load(isofoxAltSpliceJunctionCsv, altSjData);
        LOGGER.info(" Loaded {} novel splice junctions from {}", novelSpliceJunctions.size(), isofoxAltSpliceJunctionCsv);

        return ImmutableIsofoxData.builder()
                .summary(summary)
                .geneExpressions(geneExpressions)
                .fusions(fusions)
                .novelSpliceJunctions(novelSpliceJunctions)
                .build();
    }
}
