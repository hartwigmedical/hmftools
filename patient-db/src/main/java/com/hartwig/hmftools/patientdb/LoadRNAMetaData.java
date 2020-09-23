package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.SQLException;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadRNAMetaData {

    private static final Logger LOGGER = LogManager.getLogger(LoadRNAMetaData.class);

    private static final String RNA_SAMPLES_TSV = "rna_samples_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);

        String rnaSamplesTsv = cmd.getOptionValue(RNA_SAMPLES_TSV);

        if (Utils.anyNull(rnaSamplesTsv)) {
            printUsageAndExit(options);
        }

        DatabaseAccess dbAccess = databaseAccess(cmd);

        LOGGER.info("Reading RNA samples from {}", rnaSamplesTsv);
        Set<String> samples = Sets.newHashSet(Files.readAllLines(new File(rnaSamplesTsv).toPath()));
        LOGGER.info(" Loaded {} unique samples", samples.size());

        LOGGER.info("Persisting to db");
        dbAccess.writeRNA(samples);
        LOGGER.info("Complete");
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("patient-db - load rna metadata", options);
        System.exit(1);
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(RNA_SAMPLES_TSV, true, "RNA samples TSV.");
        addDatabaseCmdLineArgs(options);
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }
}
