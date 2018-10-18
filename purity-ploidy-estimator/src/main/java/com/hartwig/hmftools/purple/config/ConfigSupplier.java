package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultIntValue;
import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;
import static com.hartwig.hmftools.purple.config.StructuralVariantConfig.createStructuralVariantConfig;

import java.io.File;
import java.io.IOException;
import java.sql.Ref;
import java.util.Optional;

import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.refgenome.RefGenome;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ConfigSupplier {

    private static final Logger LOGGER = LogManager.getLogger(CommonConfig.class);

    private static final String REF_SAMPLE = "ref_sample";
    private static final String TUMOR_SAMPLE = "tumor_sample";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String OUTPUT_DIRECTORY_DEFAULT = "purple";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";
    private static final String GC_PROFILE = "gc_profile";
    private static final String REF_GENOME = "ref_genome";



    private static final String BAF = "baf";
    private static final String CIRCOS = "circos";

    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT = "min_diploid_tumor_ratio_count";
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT = 30;

    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE = "min_diploid_tumor_ratio_count_centromere";
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT = 50;

    public static void addOptions(@NotNull Options options) {
        options.addOption(REF_SAMPLE, true, "The reference sample name. Defaults to value in metadata.");
        options.addOption(TUMOR_SAMPLE, true, "The tumor sample name. Defaults to value in metadata.");
        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run.");
        options.addOption(OUTPUT_DIRECTORY, true, "The output path. Defaults to run_dir/purple/");
        options.addOption(REF_GENOME, true, "Reference genome to use. Must be one of \"hg37\" or \"hg38\".");


        options.addOption(BAF, true, "Baf file location.");
        options.addOption(CIRCOS, true, "Location of circos binary.");
        options.addOption(AMBER, true, "AMBER directory. Defaults to <run_dir>/amber");
        options.addOption(COBALT, true, "COBALT directory. Defaults to <run_dir>/cobalt");
        options.addOption(GC_PROFILE, true, "Location of GC Profile.");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT,
                true,
                "Minimum ratio count while smoothing before diploid regions become suspect.");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                true,
                "Minimum ratio count while smoothing before diploid regions become suspect while approaching centromere.");

        DBConfig.addOptions(options);
        FittingConfig.addOptions(options);
        FitScoreConfig.addOptions(options);
        SomaticConfig.addOptions(options);
        StructuralVariantConfig.addOptions(options);
    }

    private final CommonConfig commonConfig;
    private final SomaticConfig somaticConfig;
    private final StructuralVariantConfig structuralVariantConfig;
    private final CircosConfig circosConfig;
    private final DBConfig dbConfig;
    private final FittingConfig fittingConfig;
    private final SmoothingConfig smoothingConfig;
    private final FitScoreConfig fitScoreConfig;



    public ConfigSupplier(@NotNull CommandLine cmd, @NotNull Options opt) throws ParseException, IOException {
        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        if (runDirectory == null) {
            printHelp(opt);
            throw new ParseException(RUN_DIRECTORY + " is a mandatory argument");
        }

        final String gcProfile = cmd.getOptionValue(GC_PROFILE);
        if (gcProfile == null) {
            printHelp(opt);
            throw new ParseException(GC_PROFILE + " is a mandatory argument");
        }

        final String refSample;
        final String tumorSample;
        if (cmd.hasOption(REF_SAMPLE) && cmd.hasOption(TUMOR_SAMPLE)) {
            refSample = cmd.getOptionValue(REF_SAMPLE);
            tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
        } else {
            final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDirectory);
            refSample = cmd.hasOption(REF_SAMPLE) ? cmd.getOptionValue(REF_SAMPLE) : runContext.refSample();
            tumorSample = cmd.hasOption(TUMOR_SAMPLE) ? cmd.getOptionValue(TUMOR_SAMPLE) : runContext.tumorSample();
        }

        final String outputDirectory = defaultValue(cmd, OUTPUT_DIRECTORY, runDirectory + File.separator + OUTPUT_DIRECTORY_DEFAULT);
        final File outputDir = new File(outputDirectory);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + outputDirectory);
        }

        final RefGenome refGenome;
        try {
            refGenome = RefGenome.valueOf(cmd.getOptionValue(REF_GENOME).toUpperCase());
        } catch (Exception exception) {
            printHelp(opt);
            throw new ParseException(REF_GENOME + " is a mandatory argument");
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
                .refGenome(refGenome)
                .build();

        LOGGER.info("Reference Sample: {}, Tumor Sample: {}", commonConfig.refSample(), commonConfig.tumorSample());
        LOGGER.info("Output Directory: {}", commonConfig.outputDirectory());

        smoothingConfig = ImmutableSmoothingConfig.builder()
                .minDiploidTumorRatioCount(defaultIntValue(cmd, MIN_DIPLOID_TUMOR_RATIO_COUNT, MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT))
                .minDiploidTumorRatioCountAtCentromere(defaultIntValue(cmd,
                        MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                        MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT))
                .build();

        circosConfig = createCircosConfig(cmd, commonConfig);
        dbConfig = DBConfig.createConfig(cmd);
        fittingConfig = FittingConfig.createConfig(cmd);
        fitScoreConfig = FitScoreConfig.createConfig(cmd);
        somaticConfig = SomaticConfig.createSomaticConfig(cmd);
        structuralVariantConfig = createStructuralVariantConfig(cmd, opt);
    }

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
    public CircosConfig circosConfig() {
        return circosConfig;
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



    @NotNull
    private static CircosConfig createCircosConfig(@NotNull CommandLine cmd, @NotNull CommonConfig config) {
        return ImmutableCircosConfig.builder()
                .plotDirectory(config.outputDirectory() + File.separator + "plot")
                .circosDirectory(config.outputDirectory() + File.separator + "circos")
                .circosBinary(cmd.hasOption(CIRCOS) ? Optional.of(cmd.getOptionValue(CIRCOS)) : Optional.empty())
                .build();
    }

    private static void printHelp(@NotNull Options opt) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Purity Ploidy Estimator (PURPLE)", opt);
    }
}
