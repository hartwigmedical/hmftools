package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocationFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadSignatures
{
    private static final String SAMPLE_DIR = "sample_dir";

    public static void main(@NotNull String[] args) throws ParseException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        logVersion();

        try (DatabaseAccess dbAccess = createDatabaseAccess(cmd))
        {
            if(dbAccess == null)
            {
                LOGGER.error("Failed to create DB connection");
                System.exit(1);
            }

            String sampleId = cmd.getOptionValue(SAMPLE);
            String sampleFile = cmd.getOptionValue(SAMPLE_DIR);

            loadSignatureData(dbAccess, sampleId, sampleFile);

            LOGGER.info("signature allocation loading complete");
        }
        catch(Exception e)
        {
            LOGGER.error("Failed to load signature allocations", e);
            System.exit(1);
        }
    }

    private static void loadSignatureData(final DatabaseAccess dbAccess, final String sampleId, final String sampleFile)
    {
        try
        {
            final List<SignatureAllocation> sigAllocations =
                    SignatureAllocationFile.read(sampleFile);

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

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(SAMPLE_DIR, true, "Directory to read signature data from");

        return options;
    }
}
