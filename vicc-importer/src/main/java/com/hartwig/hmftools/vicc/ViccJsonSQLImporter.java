package com.hartwig.hmftools.vicc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.vicc.dao.ViccDAO;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ViccJsonSQLImporter {

    private static final Logger LOGGER = LogManager.getLogger(ViccJsonSQLImporter.class);

    private static final String VICC_JSON = "vicc_json";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String SKIP_DATABASE_WRITING = "skip_database_writing";

    private static final String VERSION = ViccJsonSQLImporter.class.getPackage().getImplementationVersion();

    public static void main(final String... args) throws ParseException, IOException, SQLException {
        LOGGER.info("Running VICC Importer v{}", VERSION);

        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        if (!validInput(cmd)) {
            printUsageAndExit(options);
        }

        String viccJsonPath = cmd.getOptionValue(VICC_JSON);
        LOGGER.info("Loading up VICC json file into memory from {}", viccJsonPath);
        List<ViccEntry> viccEntries = ViccJsonReader.buildProductionReader().readAll(viccJsonPath);
        LOGGER.info(" Loaded {} VICC entries from file.", viccEntries.size());

        if (cmd.hasOption(SKIP_DATABASE_WRITING)) {
            LOGGER.info("Skipping DB writing.");
        } else {
            ViccDAO viccDAO = connect(cmd);
            LOGGER.info("Deleting all from VICC db");
            viccDAO.deleteAll();
            LOGGER.info("Starting insertion of all VICC entries");
            int count = 0;
            for (ViccEntry viccEntry : viccEntries) {
                viccDAO.writeViccEntry(viccEntry);
                count++;
                if (count % 1000 == 0) {
                    LOGGER.info(" Completed inserting {} VICC entries into VICC db", count);
                }
            }
            LOGGER.info("Done inserting {} entries into VICC db", viccEntries.size());
        }
    }

    private static boolean validInput(@NotNull CommandLine cmd) {
        String viccJsonPath = cmd.getOptionValue(VICC_JSON);
        if (viccJsonPath == null || !pathExists(viccJsonPath)) {
            LOGGER.warn("{} has to be an existing file: {}", VICC_JSON, viccJsonPath);
            return false;
        }

        return true;
    }

    private static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    @NotNull
    private static ViccDAO connect(@NotNull CommandLine cmd) throws SQLException {
        return ViccDAO.connectToViccDAO(cmd.getOptionValue(DB_USER), cmd.getOptionValue(DB_PASS), "jdbc:" + cmd.getOptionValue(DB_URL));
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(VICC_JSON, true, "Path towards the vicc json input ");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        options.addOption(SKIP_DATABASE_WRITING, false, "If this flag is set to true, we skip the writing to the database");

        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        new HelpFormatter().printHelp("vicc-importer", options);
        System.exit(1);
    }
}
