package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV;
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion;
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;


public class LoadEvicenceData {


    private static final Logger LOGGER = LogManager.getLogger(LoadEvicenceData.class);

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    private static final String RUN_DIR = "run_dir";

    public static void main(@NotNull final String[] args) throws ParseException, SQLException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);
        final DatabaseAccess dbAccess = databaseAccess(cmd);

        final String runDirectoryPath = cmd.getOptionValue(RUN_DIR);

        if (Utils.anyNull(userName, password, databaseUrl, runDirectoryPath)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load evidence data", options);
        } else {
            final File runDirectory = new File(runDirectoryPath);
            if (runDirectory.isDirectory()) {
                final String jdbcUrl = "jdbc:" + databaseUrl;
                final DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);

                RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDirectory.toPath().toString());
                final String sample = runContext.tumorSample();

                LOGGER.info("Reading clinical Data from DB");
                final Stream<Pair<String, String>> tumorLocation = dbAccess.sampleAndTumorLocation(sample);

                LOGGER.info("Reading somatic variants from DB");
                final Stream<PotentialActionableVariant> somaticVariants = dbAccess.potentiallyActionableVariants(sample);

                LOGGER.info("Reading gene copy from DB");
                final Stream<PotentialActionableCNV> geneCopyNumber = dbAccess.potentiallyActionableCNVs(sample);

                LOGGER.info("Reading gene fusions from DB");
                final Stream<PotentialActionableFusion> geneFusion = dbAccess.potentiallyActionableFusions(sample);
            }
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(RUN_DIR, true, "Path towards the folder containing patient runs.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
