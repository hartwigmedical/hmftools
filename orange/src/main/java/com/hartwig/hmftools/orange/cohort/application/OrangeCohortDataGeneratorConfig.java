package com.hartwig.hmftools.orange.cohort.application;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.orange.util.Config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeCohortDataGeneratorConfig {

    Logger LOGGER = LogManager.getLogger(OrangeCohortDataGeneratorConfig.class);

    String DOID_JSON = "doid_json";
    String COHORT_MAPPING_TSV = "cohort_mapping_tsv";
    String OUTPUT_DIRECTORY = "output_directory";

    String LOG_DEBUG = "log_debug";

    @NotNull
    static Options createOptions() {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);

        options.addOption(DOID_JSON, true, "Path to JSON file containing the full DOID tree.");
        options.addOption(COHORT_MAPPING_TSV, true, "The cohort mapping TSV");
        options.addOption(OUTPUT_DIRECTORY, true, "Directory where output will be written to");

        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");

        return options;
    }

    @NotNull
    String doidJson();

    @NotNull
    String cohortMappingTsv();

    @NotNull
    String outputDirectory();

    @NotNull
    static OrangeCohortDataGeneratorConfig createConfig(@NotNull CommandLine cmd) throws ParseException, IOException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
            LOGGER.debug("Switched root level logging to DEBUG");
        }

        return ImmutableOrangeCohortDataGeneratorConfig.builder()
                .doidJson(nonOptionalFile(cmd, DOID_JSON))
                .cohortMappingTsv(nonOptionalFile(cmd, COHORT_MAPPING_TSV))
                .outputDirectory(outputDir(cmd, OUTPUT_DIRECTORY))
                .build();
    }

    @NotNull
    private static String nonOptionalFile(@NotNull CommandLine cmd, @NotNull String param)   {
        return Config.fileIfExists(nonOptionalValue(cmd, param));
    }

    @NotNull
    private static String nonOptionalValue(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            throw new IllegalArgumentException("Parameter must be provided: " + param);
        }

        return value;
    }

    @NotNull
    private static String outputDir(@NotNull CommandLine cmd, @NotNull String param) throws IOException {
        String value = nonOptionalValue(cmd, param);
        File outputDir = new File(value);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write to directory " + value);
        }
        return value;
    }

}
