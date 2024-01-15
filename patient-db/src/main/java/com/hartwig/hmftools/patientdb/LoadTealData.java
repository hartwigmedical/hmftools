package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
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

import com.hartwig.hmftools.common.teal.TelomereLength;
import com.hartwig.hmftools.common.teal.TelomereLengthFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LoadTealData
{
    private static final String TEAL_DIR = "teal";

    public static void main(@NotNull String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);
        configBuilder.addConfigItem(TEAL_DIR, true, "Teal output directory");
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
            String referenceId = configBuilder.getValue(REFERENCE);
            String tealDir = checkAddDirSeparator(configBuilder.getValue(TEAL_DIR));

            if(!Files.exists(Paths.get(tealDir)))
            {
                LOGGER.error("invalid Teal data directory({})", tealDir);
                System.exit(1);
            }

            LOGGER.info("loading sample({}) Teal data from {}", sampleId, tealDir);

            dbAccess.context().transaction(tr -> loadTealData(sampleId, referenceId, dbAccess, tealDir));

            LOGGER.info("Teal data loading complete");
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load Teal data: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static void loadTealData(
            String sampleId, @Nullable String referenceId, final DatabaseAccess dbAccess, final String tealDir)
    {
        @Nullable TelomereLength germlineTelomereLength = null;
        if(referenceId != null)
        {
            germlineTelomereLength = TelomereLengthFile.read(TelomereLengthFile.generateFilename(tealDir, referenceId));
        }
        TelomereLength somaticTelomereLength = TelomereLengthFile.read(TelomereLengthFile.generateFilename(tealDir, sampleId));
        dbAccess.writeTelomereLength(sampleId, germlineTelomereLength, somaticTelomereLength);
    }
}
