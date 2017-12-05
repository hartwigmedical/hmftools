package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultIntValue;
import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;

import java.io.File;
import java.util.Optional;

import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ConfigSupplier {

    private static final Logger LOGGER = LogManager.getLogger(CommonConfig.class);

    private static final String FORCE = "force";
    private static final String REF_SAMPLE = "ref_sample";
    private static final String TUMOR_SAMPLE = "tumor_sample";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String OUTPUT_DIRECTORY_DEFAULT = "purple";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";
    private static final String GC_PROFILE = "gc_profile";

    private static final String STRUCTURAL_VARIANTS = "structural_vcf";
    private static final String SOMATIC_VARIANTS = "somatic_vcf";
    private static final String SOMATIC_MIN_PEAK = "somatic_min_peak";
    private static final int SOMATIC_MIN_PEAK_DEFAULT = 50;
    private static final String SOMATIC_MIN_TOTAL = "somatic_min_total";
    private static final int SOMATIC_MIN_TOTAL_DEFAULT = 1000;

    private static final String DB_ENABLED = "db_enabled";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String BAF = "baf";
    private static final String CIRCOS = "circos";

    private static final String MIN_PURITY = "min_purity";
    private static final String MAX_PURITY = "max_purity";
    private static final String PURITY_INCREMENT = "purity_increment";
    private static final String MIN_NORM_FACTOR = "min_norm_factor";
    private static final String MAX_NORM_FACTOR = "max_norm_factor";
    private static final String NORM_FACTOR_INCREMENTS = "norm_factor_increment";

    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT = "min_diploid_tumor_ratio_count";
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT = 30;

    private static final String MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE = "min_diploid_tumor_ratio_count_centromere";
    private static final int MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT = 50;

    static final double MIN_PURITY_DEFAULT = 0.08;
    static final double MAX_PURITY_DEFAULT = 1.0;
    static final double MIN_NORM_FACTOR_DEFAULT = 0.33;
    static final double MAX_NORM_FACTOR_DEFAULT = 2.0;
    private static final double PURITY_INCREMENT_DEFAULT = 0.01;
    private static final double NORM_FACTOR_INCREMENTS_DEFAULT = 0.01;

    public static void addOptions(Options options) {
        options.addOption(REF_SAMPLE, true, "The reference sample name. Defaults to value in metadata.");
        options.addOption(TUMOR_SAMPLE, true, "The tumor sample name. Defaults to value in metadata.");
        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run.");
        options.addOption(OUTPUT_DIRECTORY, true, "The output path. Defaults to run_dir/purple/");
        options.addOption(FORCE, false, "Force recalculation of data. Do not use cached results");

        options.addOption(STRUCTURAL_VARIANTS, true, "Optional location of structural variant vcf for more accurate segmentation.");
        options.addOption(SOMATIC_VARIANTS, true, "Optional location of somatic variant vcf to assist fitting in highly-diploid samples.");

        options.addOption(BAF, true, "Baf file location.");
        options.addOption(CIRCOS, true, "Location of circos binary.");
        options.addOption(AMBER, true, "AMBER directory. Defaults to <run_dir>/amber");
        options.addOption(COBALT, true, "COBALT directory. Defaults to <run_dir>/cobalt");
        options.addOption(GC_PROFILE, true, "Location of GC Profile.");
        options.addOption(SOMATIC_MIN_PEAK, true, "Minimum number of somatic variants to consider a peak.");
        options.addOption(SOMATIC_MIN_TOTAL, true, "Minimum number of somatic variants required to assist highly diploid fits.");

        options.addOption(DB_ENABLED, false, "Persist data to DB.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        options.addOption(MIN_PURITY, true, "Minimum purity (default 0.05)");
        options.addOption(MAX_PURITY, true, "Maximum purity (default 1.0)");
        options.addOption(PURITY_INCREMENT, true, "Purity increment (default 0.01)");

        options.addOption(MIN_NORM_FACTOR, true, "Minimum norm factor (default 0.33)");
        options.addOption(MAX_NORM_FACTOR, true, "Maximum norm factor (default 2.0)");
        options.addOption(NORM_FACTOR_INCREMENTS, true, "Norm factor increments (default 0.01)");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT, true,
                "Minimum ratio count while smoothing before diploid regions become suspect.");

        options.addOption(MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE, true,
                "Minimum ratio count while smoothing before diploid regions become suspect while approaching centromere.");
    }

    private final CommonConfig commonConfig;
    private final SomaticConfig somaticConfig;
    private final StructuralVariantConfig structuralVariantConfig;
    private final BAFConfig bafConfig;
    private final CircosConfig circosConfig;
    private final DBConfig dbConfig;
    private final FittingConfig fittingConfig;
    private final SmoothingConfig smoothingConfig;

    public ConfigSupplier(@NotNull CommandLine cmd, @NotNull Options opt) throws ParseException, HartwigException {
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
        final String amberDirectory = cmd.hasOption(AMBER) ? cmd.getOptionValue(AMBER) : runDirectory + File.separator + "amber";
        final String cobaltDirectory = cmd.hasOption(COBALT) ? cmd.getOptionValue(COBALT) : runDirectory + File.separator + "cobalt";

        commonConfig = ImmutableCommonConfig.builder()
                .refSample(refSample)
                .tumorSample(tumorSample)
                .outputDirectory(outputDirectory)
                .amberDirectory(amberDirectory)
                .cobaltDirectory(cobaltDirectory)
                .forceSegmentation(cmd.hasOption(FORCE))
                .gcProfile(gcProfile)
                .build();

        LOGGER.info("Reference Sample: {}, Tumor Sample: {}", commonConfig.refSample(), commonConfig.tumorSample());
        LOGGER.info("Output Directory: {}", commonConfig.outputDirectory());

        somaticConfig = createSomaticConfig(cmd, opt);
        if (!somaticConfig.file().isPresent()) {
            LOGGER.info("No somatic vcf supplied");
        }

        structuralVariantConfig = createStructuralVariantConfig(cmd, opt);
        if (!structuralVariantConfig.file().isPresent()) {
            LOGGER.info("No structural vcf supplied");
        }

        smoothingConfig = ImmutableSmoothingConfig.builder()
                .minDiploidTumorRatioCount(defaultIntValue(cmd, MIN_DIPLOID_TUMOR_RATIO_COUNT, MIN_DIPLOID_TUMOR_RATIO_COUNT_DEFAULT))
                .minDiploidTumorRatioCountAtCentromere(defaultIntValue(cmd, MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE,
                        MIN_DIPLOID_TUMOR_RATIO_COUNT_AT_CENTROMERE_DEFAULT))
                .build();

        bafConfig = createBAFConfig(cmd, opt, commonConfig);
        circosConfig = createCircosConfig(cmd, commonConfig);
        dbConfig = createDBConfig(cmd);
        fittingConfig = createFittingConfig(cmd);
    }

    public CommonConfig commonConfig() {
        return commonConfig;
    }

    public SomaticConfig somaticConfig() {
        return somaticConfig;
    }

    public StructuralVariantConfig structuralVariantConfig() {
        return structuralVariantConfig;
    }

    public BAFConfig bafConfig() {
        return bafConfig;
    }

    public CircosConfig circosConfig() {
        return circosConfig;
    }

    public DBConfig dbConfig() {
        return dbConfig;
    }

    public FittingConfig fittingConfig() {
        return fittingConfig;
    }

    public SmoothingConfig smoothingConfig() {
        return smoothingConfig;
    }

    @NotNull
    private static SomaticConfig createSomaticConfig(@NotNull CommandLine cmd, @NotNull Options opt) throws ParseException {
        final Optional<File> file;
        if (cmd.hasOption(SOMATIC_VARIANTS)) {
            final String somaticFilename = cmd.getOptionValue(SOMATIC_VARIANTS);
            final File somaticFile = new File(somaticFilename);
            if (!somaticFile.exists()) {
                printHelp(opt);
                throw new ParseException("Unable to read somatic variants from: " + somaticFilename);
            }
            file = Optional.of(somaticFile);
        } else {
            file = Optional.empty();
        }

        return ImmutableSomaticConfig.builder()
                .file(file)
                .minTotalVariants(defaultIntValue(cmd, SOMATIC_MIN_TOTAL, SOMATIC_MIN_TOTAL_DEFAULT))
                .minPeakVariants(defaultIntValue(cmd, SOMATIC_MIN_PEAK, SOMATIC_MIN_PEAK_DEFAULT))
                .build();
    }

    @NotNull
    private static StructuralVariantConfig createStructuralVariantConfig(@NotNull CommandLine cmd, @NotNull Options opt)
            throws ParseException {
        final Optional<File> file;
        if (cmd.hasOption(STRUCTURAL_VARIANTS)) {
            final String somaticFilename = cmd.getOptionValue(STRUCTURAL_VARIANTS);
            final File somaticFile = new File(somaticFilename);
            if (!somaticFile.exists()) {
                printHelp(opt);
                throw new ParseException("Unable to read structural variants from: " + somaticFilename);
            }
            file = Optional.of(somaticFile);
        } else {
            file = Optional.empty();
        }

        return ImmutableStructuralVariantConfig.builder().file(file).build();
    }

    @NotNull
    private static CircosConfig createCircosConfig(@NotNull CommandLine cmd, @NotNull CommonConfig config) {
        return ImmutableCircosConfig.builder()
                .plotDirectory(config.outputDirectory() + File.separator + "plot")
                .circosDirectory(config.outputDirectory() + File.separator + "circos")
                .circosBinary(cmd.hasOption(CIRCOS) ? Optional.of(cmd.getOptionValue(CIRCOS)) : Optional.empty())
                .build();
    }

    @NotNull
    private static BAFConfig createBAFConfig(@NotNull final CommandLine cmd, Options opt, @NotNull final CommonConfig config)
            throws ParseException {
        if (cmd.hasOption(BAF)) {
            final String filename = cmd.getOptionValue(BAF);
            final File file = new File(filename);
            if (!file.exists()) {
                printHelp(opt);
                throw new ParseException("Unable to read bafs from: " + filename);
            }
            return ImmutableBAFConfig.builder().bafFile(file).build();
        }

        final String amberBaf = AmberBAFFile.generateAmberFilename(config.amberDirectory(), config.tumorSample());
        final File amberFile = new File(amberBaf);
        if (amberFile.exists()) {
            return ImmutableBAFConfig.builder().bafFile(amberFile).build();
        }

        final String purpleBaf = AmberBAFFile.generatePurpleFilename(config.outputDirectory(), config.tumorSample());
        final File purpleFile = new File(purpleBaf);
        if (purpleFile.exists()) {
            return ImmutableBAFConfig.builder().bafFile(purpleFile).build();
        }

        printHelp(opt);
        throw new ParseException("Baf file " + amberBaf + " not found. Please supply -baf argument.");
    }

    @NotNull
    private static DBConfig createDBConfig(@NotNull final CommandLine cmd) {
        boolean enabled = cmd.hasOption(DB_ENABLED);
        return ImmutableDBConfig.builder()
                .enabled(enabled)
                .user(enabled ? cmd.getOptionValue(DB_USER) : "")
                .password(enabled ? cmd.getOptionValue(DB_PASS) : "")
                .url("jdbc:" + (enabled ? cmd.getOptionValue(DB_URL) : ""))
                .build();
    }

    @NotNull
    private static FittingConfig createFittingConfig(@NotNull final CommandLine cmd) {
        final double minPurity = defaultValue(cmd, MIN_PURITY, MIN_PURITY_DEFAULT);
        final double maxPurity = defaultValue(cmd, MAX_PURITY, MAX_PURITY_DEFAULT);
        final double purityIncrement = defaultValue(cmd, PURITY_INCREMENT, PURITY_INCREMENT_DEFAULT);
        final double minNormFactor = defaultValue(cmd, MIN_NORM_FACTOR, MIN_NORM_FACTOR_DEFAULT);
        final double maxNormFactor = defaultValue(cmd, MAX_NORM_FACTOR, MAX_NORM_FACTOR_DEFAULT);
        final double normFactorIncrement = defaultValue(cmd, NORM_FACTOR_INCREMENTS, NORM_FACTOR_INCREMENTS_DEFAULT);

        return ImmutableFittingConfig.builder()
                .minPurity(minPurity)
                .maxPurity(maxPurity)
                .purityIncrement(purityIncrement)
                .minNormFactor(minNormFactor)
                .maxNormFactor(maxNormFactor)
                .normFactorIncrement(normFactorIncrement)
                .build();

    }

    private static void printHelp(@NotNull Options opt) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Purity Ploidy Estimator (PURPLE)", opt);
    }
}
