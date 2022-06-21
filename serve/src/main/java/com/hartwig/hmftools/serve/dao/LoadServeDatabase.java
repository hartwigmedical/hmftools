package com.hartwig.hmftools.serve.dao;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.SQLException;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import com.hartwig.hmftools.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.actionability.ActionableEventsLoader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadServeDatabase {

    private static final Logger LOGGER = LogManager.getLogger(LoadServeDatabase.class);

    private static final String SERVE_ACTIONABILITY_DIRECTORY = "serve_actionability_dir";


    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String serveActionabilityDir = nonOptionalDir(cmd, SERVE_ACTIONABILITY_DIRECTORY);
        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(RefGenomeVersion.REF_GENOME_VERSION);

        ActionableEvents actionableEvents = ActionableEventsLoader.readFromDir(serveActionabilityDir, refGenomeVersion);

        KnowledgebaseDatabaseAccess dbWriter = KnowledgebaseDatabaseAccess.databaseAccess(cmd);

        dbWriter.writeKnowledgebaseDAO(actionableEvents);
        LOGGER.info("Written serve output to database");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SERVE_ACTIONABILITY_DIRECTORY, true, "Path towards the SERVE actionability directory.");
        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, "Ref genome version to use (either '37' or '38')");

        KnowledgebaseDatabaseAccess.addDatabaseCmdLineArgs(options);
        return options;
    }

    @NotNull
    static String nonOptionalDir(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = nonOptionalValue(cmd, param);

        if (!pathExists(value) || !pathIsDirectory(value)) {
            throw new ParseException("Parameter '" + param + "' must be an existing directory: " + value);
        }

        return value;
    }

    @NotNull
    static String nonOptionalValue(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            throw new ParseException("Parameter must be provided: " + param);
        }

        return value;
    }

    static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    static boolean pathIsDirectory(@NotNull String path) {
        return Files.isDirectory(new File(path).toPath());
    }
}