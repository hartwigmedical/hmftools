package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadFlagstatData {

    private static final Logger LOGGER = LogManager.getLogger(LoadFlagstatData.class);

    private static final String SAMPLE = "sample";

    private static final String REF_FLAGSTAT_FILE = "ref_flagstat_file";
    private static final String TUMOR_FLAGSTAT_FILE = "tumor_flagstat_file";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);
        String refFlagstatFile = cmd.getOptionValue(REF_FLAGSTAT_FILE);
        String tumorFlagstatFile = cmd.getOptionValue(TUMOR_FLAGSTAT_FILE);

        if (Utils.anyNull(sample, refFlagstatFile, tumorFlagstatFile)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load Flagstat Data", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        LOGGER.info("Extracting and writing flagstats for {}", sample);

        Flagstat refFlagstat = FlagstatFile.read(refFlagstatFile);
        LOGGER.info(" Read reference sample flagstats from {}", refFlagstatFile);
        Flagstat tumorFlagstat = FlagstatFile.read(tumorFlagstatFile);
        LOGGER.info(" Read tumor sample flagstats from {}", tumorFlagstatFile);

        dbWriter.writeFlagstats(sample, refFlagstat, tumorFlagstat);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Sample for which we are going to load the flagstats");
        options.addOption(REF_FLAGSTAT_FILE, true, "Path towards the flagstat file holding the ref sample flagstats");
        options.addOption(TUMOR_FLAGSTAT_FILE, true, "Path towards the flagstat file holding the tumor sample flagstats");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
