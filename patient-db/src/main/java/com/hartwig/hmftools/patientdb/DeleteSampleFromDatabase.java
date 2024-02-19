package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.sql.SQLException;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DeleteSampleFromDatabase
{

    private static final Logger LOGGER = LogManager.getLogger(DeleteSampleFromDatabase.class);

    public static void main(@NotNull String[] args) throws ParseException, SQLException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        try(DatabaseAccess dbAccess = databaseAccess(cmd))
        {
            String sample = cmd.getOptionValue(SAMPLE);

            LOGGER.info("Removing sample '{}' from database", sample);
            dbAccess.deletePipelineDataForSample(sample);
            LOGGER.info("Complete");
        }
        catch(Exception e)
        {
            LOGGER.error("Failed to delete sample", e);
            System.exit(1);
        }
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample to delete.");
        addDatabaseCmdLineArgs(options);
        return options;
    }
}
