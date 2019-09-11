package com.hartwig.hmftools.vicc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.SQLException;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.dao.ViccDAO;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ViccJsonToSQLImporter {
    private static final Logger LOGGER = LogManager.getLogger(ViccJsonToSQLImporter.class);

    private static final String VICC_FILE = "vicc_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(final String... args) throws ParseException, IOException, SQLException {
        LOGGER.info("Attempting to load up the VICC json all file into a sql database");

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (validInput(cmd)) {
            List<ViccEntry> viccEntries = ViccJsonReader.readViccKnowledgebaseJsonFile(VICC_FILE);
            analyzeViccEntries(viccEntries);

            LOGGER.info("DONE!");
            ViccDAO viccDAO = ViccDAO.connectToViccDAO(cmd.getOptionValue(DB_USER),
                    cmd.getOptionValue(DB_PASS),
                    "jdbc:" + cmd.getOptionValue(DB_URL));
            viccDAO.deleteAll();
            int count = 0;
            for (ViccEntry viccEntry : viccEntries) {
                viccDAO.writeViccEntry(viccEntry);
                count++;
                if (count % 1000 == 0) {
                    LOGGER.info("Completed inserting " + count + " VICC entries into VICC db");
                }
            }
        } else {
            printUsageAndExit(options);
        }
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("vicc-knowledgebase-importer", options);
        System.exit(1);
    }

    private static boolean validInput(@NotNull CommandLine cmd) {
        return fileExists(cmd, VICC_FILE);
    }

    private static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value)) {
            LOGGER.warn(param + " has to be an existing file: " + value);
            return false;
        }

        return true;
    }

    private static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(VICC_FILE, true, "Path towards the vicc file ");

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

    private static void analyzeViccEntries(@NotNull List<ViccEntry> viccEntries) {
        LOGGER.info("Analyzing VICC entries - total count = " + viccEntries.size());

        LOGGER.info("Analyzing entries per SOURCE");
        countPerCategory(viccEntries, entry -> Lists.newArrayList(entry.source()));

        LOGGER.info("Analyzing entries per TAG");
        countPerCategory(viccEntries, ViccEntry::tags);

        LOGGER.info("Analyzing entries per DEVTAG");
        countPerCategory(viccEntries, ViccEntry::devTags);

        LOGGER.info("Analyzing entries per GENE");
        countPerCategory(viccEntries, ViccEntry::genes);
    }

    private static void countPerCategory(@NotNull List<ViccEntry> viccEntries, @NotNull Function<ViccEntry, List<String>> keyGenerator) {
        Map<String, Integer> countsPerKey = Maps.newHashMap();
        for (ViccEntry entry : viccEntries) {
            List<String> keys = keyGenerator.apply(entry);
            for (String key : keys) {
                countsPerKey.merge(key, 1, (a, b) -> a + b);
            }
        }

        Set<Integer> sortedCounts = Sets.newTreeSet(Comparator.reverseOrder());
        sortedCounts.addAll(countsPerKey.values());

        LOGGER.info(" Unique keys found: " + countsPerKey.keySet().size());
        for (Integer count : sortedCounts) {
            for (Map.Entry<String, Integer> entry : countsPerKey.entrySet()) {
                if (entry.getValue().equals(count)) {
                    LOGGER.info("  " + entry.getKey() + ": " + entry.getValue() + " entries");
                }
            }
        }
    }
}
