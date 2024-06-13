package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cider.Cdr3LocusSummary;
import com.hartwig.hmftools.common.cider.Cdr3LocusSummaryFile;
import com.hartwig.hmftools.common.cider.Cdr3Sequence;
import com.hartwig.hmftools.common.cider.Cdr3SequenceFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jetbrains.annotations.NotNull;

public class LoadCiderData
{
    private static final String CIDER_DIR = "cider";
    private static final String DB_SAMPLE = "db_sample";

    public static void main(@NotNull String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(DB_SAMPLE, "ID of the sample in the database (optional). Defaults to the sample ID.");
        configBuilder.addConfigItem(CIDER_DIR, true, "Cider output directory");
        addDatabaseCmdLineArgs(configBuilder, true);
        ConfigUtils.addLoggingOptions(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);
        logVersion();

        try (DatabaseAccess dbAccess = databaseAccess(configBuilder))
        {
            String sampleId = configBuilder.getValue(SAMPLE);
            String dbSampleId = configBuilder.hasValue(DB_SAMPLE) ? configBuilder.getValue(DB_SAMPLE) : sampleId;
            String ciderDir = checkAddDirSeparator(configBuilder.getValue(CIDER_DIR));

            if(!Files.exists(Paths.get(ciderDir)))
            {
                LOGGER.error("invalid Cider data directory({})", ciderDir);
                System.exit(1);
            }

            LOGGER.info("loading sample({}) dbSample({}) Cider data from {}", sampleId, dbSampleId, ciderDir);

            final String sample = sampleId;

            dbAccess.context().transaction(tr -> loadCiderData(dbSampleId, sample, dbAccess, ciderDir));

            LOGGER.info("Cider data loading complete");
        }
        catch(Exception e)
        {
            LOGGER.error("Failed to load Cider data", e);
            System.exit(1);
        }
    }

    private static void loadCiderData(final String dbSampleId, final String sampleId, final DatabaseAccess dbAccess, final String ciderDir)
    {
        // filter only partial and pass sequences
        List<Cdr3Sequence> cdr3Sequences = Cdr3SequenceFile.read(Cdr3SequenceFile.generateFilename(ciderDir, sampleId))
                .stream().filter(seq -> seq.filter().equals("PASS") || seq.filter().equals("PARTIAL"))
                .collect(Collectors.toList());

        dbAccess.writeCdr3Sequences(dbSampleId, cdr3Sequences);

        List<Cdr3LocusSummary> locusSummaries = Cdr3LocusSummaryFile.read(Cdr3LocusSummaryFile.generateFilename(ciderDir, sampleId));
        dbAccess.writeCdr3LocusSummaries(dbSampleId, locusSummaries);
    }
}
