package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jspecify.annotations.Nullable;

public class FindingApplicationConfig
{
    public final String OrangeJsonPath;
    public final @Nullable String ClinicalTranscriptsPath;
    public final String DriverGenePath;
    public final String FindingJsonPath;
    public final @Nullable Gender gender;

    private static final String ORANGE_JSON_PATH_ARG = "orange_json";
    private static final String CLINICAL_TRANSCRIPT_PATH_ARG = "clinical_transcript";
    private static final String DRIVER_GENE_PATH_ARG = "driver_gene";
    private static final String GENDER_ARG = "gender";
    private static final String FINDING_JSON_PATH_ARG = "finding_json";

    public FindingApplicationConfig(final ConfigBuilder configBuilder)
    {
        this(configBuilder.getValue(ORANGE_JSON_PATH_ARG),
                configBuilder.getValue(CLINICAL_TRANSCRIPT_PATH_ARG),
                configBuilder.getValue(DRIVER_GENE_PATH_ARG),
                configBuilder.hasValue(GENDER_ARG) ? Gender.valueOf(configBuilder.getValue(GENDER_ARG)) : null,
                configBuilder.getValue(FINDING_JSON_PATH_ARG));
    }

    public FindingApplicationConfig(
            final String orangeJsonPath,
            final @Nullable String clinicalTranscriptsPath,
            final String driverGenePath,
            final @Nullable Gender gender,
            final String findingJsonPath)
    {
        this.OrangeJsonPath = orangeJsonPath;
        this.ClinicalTranscriptsPath = clinicalTranscriptsPath;
        this.DriverGenePath = driverGenePath;
        this.FindingJsonPath = findingJsonPath;
        this.gender = gender;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(ORANGE_JSON_PATH_ARG, true, "Path to Orange JSON file");
        configBuilder.addPath(CLINICAL_TRANSCRIPT_PATH_ARG, false, "Path to clinical transcripts file");
        configBuilder.addPath(DRIVER_GENE_PATH_ARG, true, "Path to driver gene file");
        configBuilder.addConfigItem(GENDER_ARG, false, "MALE or FEMALE");
        configBuilder.addConfigItem(FINDING_JSON_PATH_ARG, true, "Path to output finding JSON file");

        addLoggingOptions(configBuilder);
    }
}
