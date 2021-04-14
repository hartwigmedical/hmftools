package com.hartwig.hmftools.ckb;

import java.io.File;
import java.nio.file.Files;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CkbImporterConfig {

    String CBK_DIR = "cbk_dir";

    String DB_USER = "db_user";
    String DB_PASS = "db_pass";
    String DB_URL = "db_url";

    String SKIP_DATABASE_WRITING = "skip_database_writing";
    String LOG_DEBUG = "log_debug";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        options.addOption(CBK_DIR, true, "Path towards the directory holding the ckb data");
        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        options.addOption(SKIP_DATABASE_WRITING, false, "If this flag is set, we skip writing to the database");

        return options;
    }

    @NotNull
    String cbkDir();

    @NotNull
    String dbUser();

    @NotNull
    String dbPass();

    @NotNull
    String dbUrl();

    boolean skipDatabaseWriting();

    @NotNull
    static CkbImporterConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        return ImmutableCkbImporterConfig.builder()
                .cbkDir(nonOptionalDir(cmd, CBK_DIR))
                .dbUser(nonOptionalValue(cmd, DB_USER))
                .dbPass(nonOptionalValue(cmd, DB_PASS))
                .dbUrl(nonOptionalValue(cmd, DB_URL))
                .skipDatabaseWriting(cmd.hasOption(SKIP_DATABASE_WRITING))
                .build();
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
