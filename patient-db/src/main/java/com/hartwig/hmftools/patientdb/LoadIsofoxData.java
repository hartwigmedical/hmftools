package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

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
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.IsofoxDAO;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadIsofoxData
{
    public static void main(@NotNull String[] args) throws ParseException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addConfig(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);
        logVersion();

        try(DatabaseAccess dbAccess = createDatabaseAccess(configBuilder))
        {
            if(dbAccess == null)
            {
                LOGGER.error("Failed to create DB connection");
                System.exit(1);
            }

            String sampleId = configBuilder.getValue(SAMPLE);
            String isofoxDir = configBuilder.getValue(ISOFOX_DIR_CFG);

            dbAccess.context().transaction(tr ->
            {
                loadIsofoxData(dbAccess, sampleId, isofoxDir);
            });

            LOGGER.info("Isofox data loading complete");
        }
        catch(Exception e)
        {
            LOGGER.error("Failed to load Isofox data", e);
            System.exit(1);
        }
    }

    private static void loadIsofoxData(final DatabaseAccess dbAccess, final String sampleId, final String isofoxDir) throws Exception
    {
        LOGGER.info("sample({}) loading Isofox data", sampleId);

        IsofoxDAO isofoxDAO = new IsofoxDAO(dbAccess.context());

        String statsFilename = RnaStatisticFile.generateFilename(isofoxDir, sampleId);

        List<String> lines = Files.readAllLines(Paths.get(statsFilename));

        RnaStatistics statistics = RnaStatisticFile.fromLines(lines);

        LOGGER.debug("sample({}) writing summary statistics to DB", sampleId);
        isofoxDAO.writeRnaStatistics(sampleId, statistics);

        String expressionFilename = GeneExpressionFile.generateFilename(isofoxDir, sampleId);

        List<GeneExpression> geneExpressions = GeneExpressionFile.read(expressionFilename);

        LOGGER.debug("sample({}) writing {} gene expression records to DB", sampleId, geneExpressions.size());
        isofoxDAO.writeGeneExpressions(sampleId, geneExpressions);

        String noveSpliceJuncFilename = NovelSpliceJunctionFile.generateFilename(isofoxDir, sampleId);
        List<NovelSpliceJunction> novelJunctions = NovelSpliceJunctionFile.read(noveSpliceJuncFilename);

        LOGGER.debug("sample({}) writing {} novel SJs records to DB", sampleId, novelJunctions.size());
        isofoxDAO.writeNovelSpliceJunctions(sampleId, novelJunctions);

        String fusionFilename = RnaFusionFile.generateFilename(isofoxDir, sampleId);

        List<RnaFusion> fusions = RnaFusionFile.read(fusionFilename);

        LOGGER.debug("sample({}) writing {} fusion records to DB", sampleId, fusions.size());
        isofoxDAO.writeRnaFusions(sampleId, fusions);
    }

    private static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
