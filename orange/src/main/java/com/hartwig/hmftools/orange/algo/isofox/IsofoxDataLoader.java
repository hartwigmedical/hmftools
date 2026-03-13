package com.hartwig.hmftools.orange.algo.isofox;

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
import com.hartwig.hmftools.orange.util.PathUtil;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class IsofoxDataLoader
{
    private static final Logger LOGGER = LogManager.getLogger(IsofoxDataLoader.class);

    public static IsofoxData load(final String tumorSampleId, final String isofoxDir) throws IOException
    {
        String summaryFile = PathUtil.mandatoryPath(RnaStatisticFile.generateFilename(isofoxDir, tumorSampleId));
        List<String> summaryLines = Files.readAllLines(Paths.get(summaryFile));
        RnaStatistics summary = RnaStatisticFile.fromLines(summaryLines);

        LOGGER.debug((" loaded summary from " + summaryFile));

        String isofoxGeneFile = PathUtil.mandatoryPath(GeneExpressionFile.generateFilename(isofoxDir, tumorSampleId));
        List<GeneExpression> geneExpressions = GeneExpressionFile.read(isofoxGeneFile);
        LOGGER.debug(" loaded {} gene expressions from {}", geneExpressions.size(), isofoxGeneFile);

        String fusionFile = PathUtil.mandatoryPath(RnaFusionFile.generateFilename(isofoxDir, tumorSampleId));
        List<RnaFusion> fusions = RnaFusionFile.read(fusionFile);
        LOGGER.debug(" loaded {} fusions from {}", fusions.size(), fusionFile);

        String altSpliceJunctionFile = PathUtil.mandatoryPath(NovelSpliceJunctionFile.generateFilename(isofoxDir, tumorSampleId));
        List<NovelSpliceJunction> novelSpliceJunctions = NovelSpliceJunctionFile.read(altSpliceJunctionFile);
        LOGGER.debug(" loaded {} novel splice junctions from {}", novelSpliceJunctions.size(), altSpliceJunctionFile);

        return ImmutableIsofoxData.builder()
                .summary(summary)
                .geneExpressions(geneExpressions)
                .fusions(fusions)
                .novelSpliceJunctions(novelSpliceJunctions)
                .build();
    }
}
