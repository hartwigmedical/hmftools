package com.hartwig.hmftools.patientdb;

import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.amber.AmberPatient;
import com.hartwig.hmftools.common.amber.AmberPatientFactory;
import com.hartwig.hmftools.common.amber.AmberSample;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RefreshAmberPatient {

    private static final Logger LOGGER = LogManager.getLogger(RefreshAmberPatient.class);

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull String[] args) throws ParseException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);

        try (final DatabaseAccess dbAccess = databaseAccess(cmd)) {

            LOGGER.info("Truncating AMBER patients");
            dbAccess.truncateAmberPatients();

            final List<AmberSample> allSamples = dbAccess.readAmberSamples();
            for (int i = 0; i < allSamples.size(); i++) {
                final AmberSample sample = allSamples.get(i);
                LOGGER.info("Processing " + (i + 1) + " of " + allSamples.size() + ": " + sample.sampleId());
                final List<AmberPatient> allPatients = dbAccess.readAmberPatients();
                final List<AmberPatient> updatedPatients = AmberPatientFactory.create(0.9, sample, allSamples, allPatients);
                dbAccess.writeAmberPatients(sample.sampleId(), updatedPatients);
            }

        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull CommandLine cmd) throws SQLException {
        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
