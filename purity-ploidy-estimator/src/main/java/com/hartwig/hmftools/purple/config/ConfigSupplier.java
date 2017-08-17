package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;

import java.io.File;
import java.util.Optional;

import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.copynumber.freec.FreecFileLoader;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.baf.TumorBAFFile;

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
    private static final String FREEC_DIRECTORY = "freec_dir";
    private static final String REF_SAMPLE = "ref_sample";
    private static final String TUMOR_SAMPLE = "tumor_sample";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String OUTPUT_DIRECTORY_DEFAULT = "purple";

    private static final String STRUCTURAL_VARIANTS = "structural_vcf";
    private static final String SOMATIC_VARIANTS = "somatic_vcf";
    private static final String BAF_VARIANTS = "baf_vcf";
    private static final String BAF = "baf";

    public static void addOptions(Options options) {
        options.addOption(REF_SAMPLE, true, "The reference sample name. Defaults to value in metadata.");
        options.addOption(TUMOR_SAMPLE, true, "The tumor sample name. Defaults to value in metadata.");
        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run.");
        options.addOption(OUTPUT_DIRECTORY, true, "The output path. Defaults to run_dir/purple/");
        options.addOption(FREEC_DIRECTORY, true, "The freec data path. Defaults to run_dir/copyNumber/refSample_tumorSample/freec/");
        options.addOption(FORCE, false, "Force recalculation of data. Do not use cached results");

        options.addOption(STRUCTURAL_VARIANTS, true, "Optional location of structural variant vcf for more accurate segmentation.");
        options.addOption(SOMATIC_VARIANTS, true, "Optional location of somatic variant vcf to assist fitting in highly-diploid samples.");

        options.addOption(BAF_VARIANTS, true, "Location of vcf to calculate BAF.");
        options.addOption(BAF, true, "Baf file location.");
    }

    private final CommonConfig commonConfig;
    private final SomaticConfig somaticConfig;
    private final StructuralVariantConfig structuralVariantConfig;
    private final BAFConfig bafConfig;

    public ConfigSupplier(CommandLine cmd, Options opt) throws ParseException, HartwigException {
        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        if (runDirectory == null) {
            printHelp(opt);
            throw new ParseException(RUN_DIRECTORY + " is a mandatory argument");
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
        final String freecDirectory = freecDirectory(cmd, runDirectory, refSample, tumorSample);
        commonConfig = new CommonConfig(refSample, tumorSample, outputDirectory, runDirectory, freecDirectory, cmd.hasOption(FORCE));

        LOGGER.info("Reference Sample: {}, Tumor Sample: {}", commonConfig.refSample(), commonConfig.tumorSample());
        LOGGER.info("Run Directory: {}", commonConfig.runDirectory());
        LOGGER.info("Output Directory: {}", commonConfig.outputDirectory());

        somaticConfig = createSomaticConfig(cmd, opt);
        if (!somaticConfig.file().isPresent()) {
            LOGGER.info("No somatic vcf supplied");
        }

        structuralVariantConfig = createStructuralVariantConfig(cmd, opt);
        if (!structuralVariantConfig.file().isPresent()) {
            LOGGER.info("No structural vcf supplied");
        }

        bafConfig = createBAFConfig(cmd, opt, commonConfig);
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

    @NotNull
    private static String freecDirectory(@NotNull final CommandLine cmd, @NotNull final String runDirectory,
            @NotNull final String refSample, @NotNull final String tumorSample) {
        return cmd.hasOption(FREEC_DIRECTORY)
                ? cmd.getOptionValue(FREEC_DIRECTORY)
                : FreecFileLoader.getFreecBasePath(runDirectory, refSample, tumorSample);
    }

    @NotNull
    private static SomaticConfig createSomaticConfig(CommandLine cmd, Options opt) throws ParseException {

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

        return ImmutableSomaticConfig.builder().file(file).build();
    }

    @NotNull
    private static StructuralVariantConfig createStructuralVariantConfig(CommandLine cmd, Options opt) throws ParseException {

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
    private static BAFConfig createBAFConfig(CommandLine cmd, Options opt, CommonConfig config) throws ParseException {

        final ImmutableBAFConfig.Builder builder = ImmutableBAFConfig.builder().bafFile(Optional.empty()).bafVCFFile(Optional.empty());

        if (cmd.hasOption(BAF)) {
            final String filename = cmd.getOptionValue(BAF);
            final File file = new File(filename);
            if (!file.exists()) {
                printHelp(opt);
                throw new ParseException("Unable to read bafs from: " + filename);
            }
            return builder.bafFile(Optional.of(file)).build();
        }

        final String cachedBafFilename = TumorBAFFile.generateFilename(config.outputDirectory(), config.tumorSample());
        final File cachedFile = new File(cachedBafFilename);
        if (cachedFile.exists()) {
            return builder.bafFile(Optional.of(cachedFile)).build();
        }

        if (cmd.hasOption(BAF_VARIANTS)) {
            final String filename = cmd.getOptionValue(BAF_VARIANTS);
            final File file = new File(filename);
            if (!file.exists()) {
                printHelp(opt);
                throw new ParseException("Unable to read bafs from: " + filename);
            }
            return builder.bafVCFFile(Optional.of(file)).build();
        }

        printHelp(opt);
        throw new ParseException("Cached baf file " + cachedBafFilename + " not found. Please supply one of -baf or -baf_vcf arguments.");
    }

    private static void printHelp(Options opt) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Purity Ploidy Estimator (PURPLE)", opt);
    }

}
