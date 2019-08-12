package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
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

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String RUN_DIR = "run_dir";
    private static final String TUMOR_SAMPLE = "tumor_sample";

    public static void main(@NotNull final String[] args) throws ParseException, SQLException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);

        final String runDirectoryPath = cmd.getOptionValue(RUN_DIR);
        final String tumor_sample = cmd.getOptionValue(TUMOR_SAMPLE);

        if (Utils.anyNull(userName, password, databaseUrl, runDirectoryPath, tumor_sample)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load chord data", options);
        } else {
            final File runDirectory = new File(runDirectoryPath);
            if (runDirectory.isDirectory()) {
                final String jdbcUrl = "jdbc:" + databaseUrl;
                final DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);

                LOGGER.info(String.format("Extracting and writing chord for %s", runDirectoryPath));
                try {
                    ChordAnalysis chordAnalysis = generateChordForRun(runDirectoryPath, tumor_sample);
                    dbWriter.writeChord(tumor_sample, chordAnalysis);

                } catch (IOException e) {
                    LOGGER.warn(String.format("Cannot extract chord for %s.", runDirectoryPath));
                }
            } else {
                if (!runDirectory.exists()) {
                    LOGGER.warn("dir " + runDirectory + " does not exist.");
                }
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("patient-db - load chord data", options);
            }
        }
    }

    @NotNull
    private static ChordAnalysis generateChordForRun(@NotNull String runDirectoryPath, @NotNull String tumorSample) throws IOException {
        String chordFile = ChordFileReader.generateFilename(runDirectoryPath, tumorSample);
        return ChordFileReader.read(chordFile);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(RUN_DIR, true, "Path towards the folder containing patient runs.");
        options.addOption(TUMOR_SAMPLE, true, "The tumor sample.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}