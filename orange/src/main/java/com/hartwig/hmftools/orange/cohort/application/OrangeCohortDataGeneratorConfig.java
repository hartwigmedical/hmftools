package com.hartwig.hmftools.orange.cohort.application;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeCohortDataGeneratorConfig
{
    String DOID_JSON = "doid_json";
    String COHORT_MAPPING_TSV = "cohort_mapping_tsv";

    static void registerConfig(final ConfigBuilder configBuilder)
    {
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(DOID_JSON, true, "Path to JSON file containing the full DOID tree.");
        configBuilder.addPath(COHORT_MAPPING_TSV, true, "The cohort mapping TSV");
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }

    @NotNull
    String doidJson();

    @NotNull
    String cohortMappingTsv();

    @NotNull
    String outputDirectory();

    @NotNull
    static OrangeCohortDataGeneratorConfig createConfig(final ConfigBuilder configBuilder) throws IOException
    {
        String outputDir = parseOutputDir(configBuilder);

        File outputDirFile = new File(outputDir);
        if(!outputDirFile.exists() && !outputDirFile.mkdirs())
        {
            throw new IOException("Unable to write to directory " + outputDir);
        }

        return ImmutableOrangeCohortDataGeneratorConfig.builder()
                .doidJson(configBuilder.getValue(DOID_JSON))
                .cohortMappingTsv(configBuilder.getValue(COHORT_MAPPING_TSV))
                .outputDirectory(outputDir)
                .build();
    }
}
