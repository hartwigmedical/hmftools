package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultIntValue;
import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;
import static com.hartwig.hmftools.purple.config.StructuralVariantConfig.createStructuralVariantConfig;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.io.exception.MalformedFileException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ConfigSupplier {

    private static final Logger LOGGER = LogManager.getLogger(CommonConfig.class);

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String OUTPUT_DIRECTORY_DEFAULT = "purple";
    private static final String GC_PROFILE = "gc_profile";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";

    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT = "min_diploid_tumor_ratio_count";
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT = 30;

    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE = "min_diploid_tumor_ratio_count_centromere";
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT = 50;

    public static void addOptions(@NotNull Options options) {
        options.addOption(REF_SAMPLE, true, "Name of the reference sample. This should correspond to the value used in AMBER and COBALT.");
        options.addOption(TUMOR_SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in AMBER and COBALT.");
        options.addOption(RUN_DIRECTORY,
                true,
                "If provided, default values of <run_dir>/amber, <run_dir>/cobalt and <run_dir>/purple will be supplied for amber, cobalt and output_dir parameters respectively.");




        options.addOption(GC_PROFILE, true, "Path to GC profile.");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT,
                true,
                "Minimum ratio count while smoothing before diploid regions become suspect.");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                true,
                "Minimum ratio count while smoothing before diploid regions become suspect while approaching centromere.");

        options.addOption(OUTPUT_DIRECTORY, true, "Path to the output directory. Required if <run_dir> not set, otherwise defaults to run_dir/purple/");
        options.addOption(COBALT, true, "Path to COBALT output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/cobalt.");
        options.addOption(AMBER, true, "Path to AMBER output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/amber");

        DBConfig.addOptions(options);
        FittingConfig.addOptions(options);
        FitScoreConfig.addOptions(options);
        SomaticConfig.addOptions(options);
        StructuralVariantConfig.addOptions(options);
        RefGenomeData.addOptions(options);
        ChartConfig.addOptions(options);
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

    private final CobaltData cobaltData;
    private final AmberData amberData;

    public ConfigSupplier(@NotNull CommandLine cmd, @NotNull Options opt) throws ParseException, IOException {

        final String gcProfile = cmd.getOptionValue(GC_PROFILE);
        if (gcProfile == null) {
            printHelp(opt);
            throw new ParseException(GC_PROFILE + " is a mandatory argument");
        }

        if (!cmd.hasOption(RUN_DIRECTORY)) {
            if (!cmd.hasOption(REF_SAMPLE)) {
                printHelp(opt);
                throw new ParseException(REF_SAMPLE + " is a mandatory argument");
            }
            if (!cmd.hasOption(TUMOR_SAMPLE)) {
                printHelp(opt);
                throw new ParseException(TUMOR_SAMPLE + " is a mandatory argument");
            }

            if (!cmd.hasOption(OUTPUT_DIRECTORY)) {
                printHelp(opt);
                throw new ParseException("Either " + RUN_DIRECTORY + " or " + OUTPUT_DIRECTORY + " must be supplied as a parameter");
            }

            if (!cmd.hasOption(AMBER)) {
                printHelp(opt);
                throw new ParseException("Either " + RUN_DIRECTORY + " or " + AMBER + " must be supplied as a parameter");
            }

            if (!cmd.hasOption(COBALT)) {
                printHelp(opt);
                throw new ParseException("Either " + RUN_DIRECTORY + " or " + COBALT + " must be supplied as a parameter");
            }
        }

        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        final String refSample;
        final String tumorSample;
        if (cmd.hasOption(REF_SAMPLE) && cmd.hasOption(TUMOR_SAMPLE)) {
            refSample = cmd.getOptionValue(REF_SAMPLE);
            tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
        } else {

            try {
                final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDirectory);
                refSample = cmd.hasOption(REF_SAMPLE) ? cmd.getOptionValue(REF_SAMPLE) : runContext.refSample();
                tumorSample = cmd.hasOption(TUMOR_SAMPLE) ? cmd.getOptionValue(TUMOR_SAMPLE) : runContext.tumorSample();
            } catch (MalformedFileException e) {
                printHelp(opt);
                throw new ParseException("Unable to resolve " + REF_SAMPLE + " and " + TUMOR_SAMPLE + " from meta data. These parameters must be supplied.");
            }
        }

        final String outputDirectory = defaultValue(cmd, OUTPUT_DIRECTORY, runDirectory + File.separator + OUTPUT_DIRECTORY_DEFAULT);
        final File outputDir = new File(outputDirectory);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + outputDirectory);
        }

        final String amberDirectory = cmd.hasOption(AMBER) ? cmd.getOptionValue(AMBER) : runDirectory + File.separator + "amber";
        final String cobaltDirectory = cmd.hasOption(COBALT) ? cmd.getOptionValue(COBALT) : runDirectory + File.separator + "cobalt";

        commonConfig = ImmutableCommonConfig.builder()
                .refSample(refSample)
                .tumorSample(tumorSample)
                .outputDirectory(outputDirectory)
                .amberDirectory(amberDirectory)
                .cobaltDirectory(cobaltDirectory)
                .gcProfile(gcProfile)
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
        somaticConfig = SomaticConfig.createSomaticConfig(cmd);
        structuralVariantConfig = createStructuralVariantConfig(cmd, opt);

        cobaltData = CobaltData.createCobaltData(commonConfig);
        amberData = AmberData.createAmberData(commonConfig);

        refGenomeData = RefGenomeData.createRefGenomeConfig(cmd, cobaltData);
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
    public ChartConfig circosConfig() {
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


    private static void printHelp(@NotNull Options opt) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Purity Ploidy Estimator (PURPLE)", opt);
    }
}
