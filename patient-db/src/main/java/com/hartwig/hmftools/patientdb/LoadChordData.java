package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LoadChordData {

    private static final Logger LOGGER = LogManager.getLogger(LoadChordData.class);

    private static final String SAMPLE = "sample";

    private static final String PREDICTION_FILE = "prediction_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);

        final String predictionFile = cmd.getOptionValue(PREDICTION_FILE);
        final String sample = cmd.getOptionValue(SAMPLE);

        if (Utils.anyNull(userName, password, databaseUrl, predictionFile, sample)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load chord data", options);
        } else {
            final File fileObject = new File(predictionFile);
            if (fileObject.isFile()) {
                final String jdbcUrl = "jdbc:" + databaseUrl;
                final DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);

                LOGGER.info(String.format("Extracting and writing chord for %s", predictionFile));
                try {
                    ChordAnalysis chordAnalysis = generateChordForRun(predictionFile);
                    dbWriter.writeChord(sample, chordAnalysis);

                } catch (IOException e) {
                    LOGGER.warn(String.format("Cannot extract chord for %s.", predictionFile));
                }
            } else {
                if (!fileObject.exists()) {
                    LOGGER.warn("file " + predictionFile + " does not exist.");
                }
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("patient-db - load chord data", options);
            }
        }
    }

    @NotNull
    private static ChordAnalysis generateChordForRun(@NotNull String chordFile) throws IOException {
        return ChordFileReader.read(chordFile);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "The tumor sample.");
        options.addOption(PREDICTION_FILE, true, "Path towards the chord prediction file.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}