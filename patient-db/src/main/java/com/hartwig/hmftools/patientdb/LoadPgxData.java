package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadPgxData {
    private static final Logger LOGGER = LogManager.getLogger(LoadPgxData.class);

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String SAMPLE = "sample";


    public static void main(@NotNull final String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);

        String sample = cmd.getOptionValue(SAMPLE);


        if (Utils.anyNull(userName, password, databaseUrl )) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load metrics data", options);
        } else {
            final String jdbcUrl = "jdbc:" + databaseUrl;
            final DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);

            LOGGER.info("Extracting and writing pgx for {}", sample);


        }

    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Sample for which we are going to load the metrics");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }


}
