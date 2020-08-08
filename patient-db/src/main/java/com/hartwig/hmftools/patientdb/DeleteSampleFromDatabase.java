package com.hartwig.hmftools.patientdb;

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

public class DeleteSampleFromDatabase {

    private static final Logger LOGGER = LogManager.getLogger(DeleteSampleFromDatabase.class);

    private static final String SAMPLE = "sample";

    public static void main(@NotNull String[] args) throws ParseException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);
        DatabaseAccess dbAccess = databaseAccess(cmd);
        String sample = cmd.getOptionValue(SAMPLE);

        LOGGER.info("Removing sample: " + sample + " from hmf database");
        dbAccess.deleteAllDataForSample(sample);
        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample to delete.");
        addDatabaseCmdLineArgs(options);
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }
}
