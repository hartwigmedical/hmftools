package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultIntValue;
import static com.hartwig.hmftools.purple.config.StructuralVariantConfig.createStructuralVariantConfig;

import java.io.File;
import java.io.IOException;
import java.util.StringJoiner;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ConfigSupplier {

    private static final Logger LOGGER = LogManager.getLogger(CommonConfig.class);

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String GC_PROFILE = "gc_profile";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";
    private static final String TUMOR_ONLY = "tumor_only";

    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT = "min_diploid_tumor_ratio_count";
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT = 30;

    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE = "min_diploid_tumor_ratio_count_centromere";
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT = 150;

    public static void addOptions(@NotNull Options options) {
        options.addOption(TUMOR_ONLY, false, "Tumor only mode. Disables somatic fitting.");
        options.addOption(REF_SAMPLE, true, "Name of the reference sample. This should correspond to the value used in AMBER and COBALT.");
        options.addOption(TUMOR_SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in AMBER and COBALT.");

        options.addOption(GC_PROFILE, true, "Path to GC profile.");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT,
                true,
                "Minimum ratio count while smoothing before diploid regions become suspect.");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                true,
                "Minimum ratio count while smoothing before diploid regions become suspect while approaching centromere.");

        options.addOption(OUTPUT_DIRECTORY,
                true,
                "Path to the output directory. Required if <run_dir> not set, otherwise defaults to run_dir/purple/");
        options.addOption(COBALT,
                true,
                "Path to COBALT output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/cobalt.");
        options.addOption(AMBER,
                true,
                "Path to AMBER output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/amber");

        DBConfig.addOptions(options);
        FittingConfig.addOptions(options);
        FitScoreConfig.addOptions(options);
        SomaticConfig.addOptions(options);
        StructuralVariantConfig.addOptions(options);
        RefGenomeData.addOptions(options);
        ChartConfig.addOptions(options);
        DriverCatalogConfig.addOptions(options);
    }

    private final CommonConfig commonConfig;
    private final SomaticConfig somaticConfig;
    private final StructuralVariantConfig structuralVariantConfig;
    private final ChartConfig chartConfig;
    private final DBConfig dbConfig;
    private final FittingConfig fittingConfig;
    private final SmoothingConfig smoothingConfig;
    private final FitScoreConfig fitScoreConfig;
    private final RefGenomeData refGenomeData;
    private final DriverCatalogConfig driverCatalogConfig;

    private final CobaltData cobaltData;
    private final AmberData amberData;

    public ConfigSupplier(@NotNull final String version, @NotNull CommandLine cmd, @NotNull Options opt)
            throws ParseException, IOException {

        final StringJoiner missingJoiner = new StringJoiner(", ");
        final String gcProfile = parameter(cmd, GC_PROFILE, missingJoiner);
        final String refSample = cmd.getOptionValue(REF_SAMPLE, "");
        final String tumorSample = parameter(cmd, TUMOR_SAMPLE, missingJoiner);
        final String outputDirectory = parameter(cmd, OUTPUT_DIRECTORY, missingJoiner);
        final String amberDirectory = parameter(cmd, AMBER, missingJoiner);
        final String cobaltDirectory = parameter(cmd, COBALT, missingJoiner);
        final String missing = missingJoiner.toString();

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        final File outputDir = new File(outputDirectory);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + outputDirectory);
        }

        commonConfig = ImmutableCommonConfig.builder()
                .version(version)
                .refSample(refSample)
                .tumorSample(tumorSample)
                .outputDirectory(outputDirectory)
                .amberDirectory(amberDirectory)
                .cobaltDirectory(cobaltDirectory)
                .gcProfile(gcProfile)
                .tumorOnly(cmd.hasOption(TUMOR_ONLY))
                .build();

        LOGGER.info("Reference Sample: {}, Tumor Sample: {}", commonConfig.refSample(), commonConfig.tumorSample());
        LOGGER.info("Output Directory: {}", commonConfig.outputDirectory());

        smoothingConfig = ImmutableSmoothingConfig.builder()
                .minDiploidTumorRatioCount(defaultIntValue(cmd, MIN_DIPLOID_TUMOR_RATIO_COUNT, MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT))
                .minDiploidTumorRatioCountAtCentromere(defaultIntValue(cmd,
                        MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                        MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT))
                .build();

        chartConfig = ChartConfig.createCircosConfig(cmd, commonConfig);
        dbConfig = DBConfig.createConfig(cmd);
        fittingConfig = FittingConfig.createConfig(cmd);
        fitScoreConfig = FitScoreConfig.createConfig(cmd);
        structuralVariantConfig = createStructuralVariantConfig(cmd, opt);
        refGenomeData = RefGenomeData.createRefGenomeConfig(cmd);
        driverCatalogConfig = DriverCatalogConfig.createConfig(cmd, refGenomeData);

        cobaltData = CobaltData.createCobaltData(commonConfig);
        amberData = AmberData.createAmberData(commonConfig);

        somaticConfig = SomaticConfig.createSomaticConfig(cmd, amberData);
    }

    @NotNull
    public CobaltData cobaltData() {
        return cobaltData;
    }

    @NotNull
    public AmberData amberData() {
        return amberData;
    }

    @NotNull
    public RefGenomeData refGenomeConfig() {
        return refGenomeData;
    }

    @NotNull
    public FitScoreConfig fitScoreConfig() {
        return fitScoreConfig;
    }

    @NotNull
    public CommonConfig commonConfig() {
        return commonConfig;
    }

    @NotNull
    public SomaticConfig somaticConfig() {
        return somaticConfig;
    }

    @NotNull
    public StructuralVariantConfig structuralVariantConfig() {
        return structuralVariantConfig;
    }

    @NotNull
    public ChartConfig chartConfig() {
        return chartConfig;
    }

    @NotNull
    public DBConfig dbConfig() {
        return dbConfig;
    }

    @NotNull
    public FittingConfig fittingConfig() {
        return fittingConfig;
    }

    @NotNull
    public SmoothingConfig smoothingConfig() {
        return smoothingConfig;
    }

    public DriverCatalogConfig driverCatalogConfig() {
        return driverCatalogConfig;
    }

    @NotNull
    static String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing) {
        final String value = cmd.getOptionValue(parameter);
        if (value == null) {
            missing.add(parameter);
            return "";
        }
        return value;
    }
}
