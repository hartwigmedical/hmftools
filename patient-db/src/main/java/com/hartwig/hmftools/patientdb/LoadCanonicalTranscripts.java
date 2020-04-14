package com.hartwig.hmftools.patientdb;

import java.sql.SQLException;

import com.hartwig.hmftools.common.genome.region.CanonicalTranscriptFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadCanonicalTranscripts {

    private static final Logger LOGGER = LogManager.getLogger(LoadCanonicalTranscripts.class);

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull String[] args) throws ParseException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        LOGGER.info("Persisting transcripts to database");
        dbAccess.writeCanonicalTranscripts("GRCh37", CanonicalTranscriptFactory.create37());
        dbAccess.writeCanonicalTranscripts("GRCh38", CanonicalTranscriptFactory.create38());

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

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
