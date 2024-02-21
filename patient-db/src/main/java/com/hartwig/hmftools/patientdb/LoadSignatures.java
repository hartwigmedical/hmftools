package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadSignatures
{
    private static final String SAMPLE_DIR = "sample_dir";

    public static void main(@NotNull String[] args) throws ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(SAMPLE_DIR, true, "Directory to read signature data from");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);
        String sampleDir = configBuilder.getValue(SAMPLE_DIR);

        try(DatabaseAccess dbAccess = createDatabaseAccess(configBuilder))
        {
            if(dbAccess == null)
            {
                LOGGER.error("Failed to create DB connection");
                System.exit(1);
            }

            loadSignatureData(dbAccess, sample, sampleDir);

            LOGGER.info("signature allocation loading complete");
        }
        catch(Exception e)
        {
            LOGGER.error("Failed to load signature allocations", e);
            System.exit(1);
        }
    }

    private static void loadSignatureData(final DatabaseAccess dbAccess, final String sampleId, final String sampleDir)
    {
        try
        {
            final List<SignatureAllocation> sigAllocations =
                    SignatureAllocationFile.read(SignatureAllocationFile.generateFilename(sampleDir, sampleId));

            if(!sigAllocations.isEmpty())
            {
                LOGGER.info("sample({}) writing {} allocations to database", sampleId, sigAllocations.size());
                dbAccess.writeSignatures(sampleId, sigAllocations);
            }
            else
            {
                LOGGER.info("sample({}) has not signature allocations", sampleId);
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load sample({}) allocations: {}", sampleId, e.toString());
        }
    }
}
