package com.hartwig.hmftools.orange.cohort.application;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.IOException;

import com.hartwig.hmftools.orange.util.Config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeCohortGeneratorConfig {

    String DOID_JSON = "doid_json";
    String COHORT_MAPPING_TSV = "cohort_mapping_tsv";
    String OUTPUT_DIRECTORY = "output_directory";

    @NotNull
    static Options createOptions() {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);

        options.addOption(DOID_JSON, true, "Path to JSON file containing the full DOID tree.");
        options.addOption(COHORT_MAPPING_TSV, true, "The cohort mapping TSV");
        options.addOption(OUTPUT_DIRECTORY, true, "Directory where output will be written to");

        return options;
    }

    @NotNull
    String doidJson();

    @NotNull
    String cohortMappingTsv();

    @NotNull
    String outputDirectory();

    @NotNull
    static OrangeCohortGeneratorConfig createConfig(@NotNull CommandLine cmd) throws ParseException, IOException {
        return ImmutableOrangeCohortGeneratorConfig.builder()
                .doidJson(Config.nonOptionalFile(cmd, DOID_JSON))
                .cohortMappingTsv(Config.nonOptionalFile(cmd, COHORT_MAPPING_TSV))
                .outputDirectory(Config.outputDir(cmd, OUTPUT_DIRECTORY))
                .build();
    }
}
