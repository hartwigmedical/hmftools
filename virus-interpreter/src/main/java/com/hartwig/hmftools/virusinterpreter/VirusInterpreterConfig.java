package com.hartwig.hmftools.virusinterpreter;

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
public interface VirusInterpreterConfig {

    String SAMPLE_ID = "sample_id";
    String VIRUS_BREAKEND_TSV = "virus_breakend_tsv";

    String TAXONOMY_DB_TSV = "taxonomy_db_tsv";
    String VIRUS_INTERPRETATION_TSV = "virus_interpretation_tsv";
    String VIRUS_BLACKLIST_TSV = "virus_blacklist_tsv";

    String OUTPUT_DIRECTORY = "output_dir";

    String LOG_DEBUG = "log_debug";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE_ID, true, "The sample ID for which virus interpreter will run.");
        options.addOption(VIRUS_BREAKEND_TSV, true, "Path towards the virus breakend TSV.");
        options.addOption(TAXONOMY_DB_TSV, true, "Path towards a TSV containing a mapping from taxid to taxonomy name.");
        options.addOption(VIRUS_INTERPRETATION_TSV, true, "Path towards a TSV containing interpretation rules for viruses.");
        options.addOption(VIRUS_BLACKLIST_TSV, true, "Path towards a TSV containing blacklisting for specific viruses.");

        options.addOption(OUTPUT_DIRECTORY, true, "Path to where the virus interpretation output will be written to.");

        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");

        return options;
    }

    @NotNull
    String sampleId();

    @NotNull
    String virusBreakendTsv();

    @NotNull
    String taxonomyDbTsv();

    @NotNull
    String virusInterpretationTsv();

    @NotNull
    String virusBlacklistTsv();

    @NotNull
    String outputDir();

    @NotNull
    static VirusInterpreterConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        return ImmutableVirusInterpreterConfig.builder().sampleId(nonOptionalValue(cmd, SAMPLE_ID))
                .virusBreakendTsv(nonOptionalFile(cmd, VIRUS_BREAKEND_TSV))
                .taxonomyDbTsv(nonOptionalFile(cmd, TAXONOMY_DB_TSV))
                .virusInterpretationTsv(nonOptionalFile(cmd, VIRUS_INTERPRETATION_TSV))
                .virusBlacklistTsv(nonOptionalFile(cmd, VIRUS_BLACKLIST_TSV))
                .outputDir(nonOptionalDir(cmd, OUTPUT_DIRECTORY))
                .build();
    }

    @NotNull
    static String nonOptionalValue(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            throw new ParseException("Parameter must be provided: " + param);
        }

        return value;
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
    static String nonOptionalFile(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = nonOptionalValue(cmd, param);

        if (!pathExists(value)) {
            throw new ParseException("Parameter '" + param + "' must be an existing file: " + value);
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
