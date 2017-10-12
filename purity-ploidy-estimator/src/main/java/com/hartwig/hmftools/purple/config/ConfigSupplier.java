package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;

import java.io.File;
import java.util.Optional;

import com.hartwig.hmftools.common.baf.TumorBAFFile;
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
    private static final String FREEC_DIRECTORY = "freec_dir";
    private static final String REF_SAMPLE = "ref_sample";
    private static final String TUMOR_SAMPLE = "tumor_sample";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String OUTPUT_DIRECTORY_DEFAULT = "purple";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";

    private static final String STRUCTURAL_VARIANTS = "structural_vcf";
    private static final String SOMATIC_VARIANTS = "somatic_vcf";
    private static final String BAF_VARIANTS = "baf_vcf";
    private static final String BAF = "baf";
    private static final String CIRCOS = "circos";

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
        options.addOption(CIRCOS, true, "Location of circos binary.");
        options.addOption(AMBER, true, "AMBER directory. Defaults to <run_dir>/amber");
        options.addOption(COBALT, true, "COBALT directory. Defaults to <run_dir>/cobalt");
    }

    private final CommonConfig commonConfig;
    private final SomaticConfig somaticConfig;
    private final StructuralVariantConfig structuralVariantConfig;
    private final BAFConfig bafConfig;
    private final CircosConfig circosConfig;

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
        final String amberDirectory = cmd.hasOption(AMBER) ? cmd.getOptionValue(AMBER) : runDirectory + File.separator + "amber";
        final String cobaltDirectory = cmd.hasOption(COBALT) ? cmd.getOptionValue(COBALT) : runDirectory + File.separator + "cobalt";

        commonConfig = ImmutableCommonConfig.builder()
                .refSample(refSample)
                .tumorSample(tumorSample)
                .outputDirectory(outputDirectory)
                .amberDirectory(amberDirectory)
                .cobaltDirectory(cobaltDirectory)
                .forceSegmentation(cmd.hasOption(FORCE))
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

        bafConfig = createBAFConfig(cmd, opt, commonConfig);
        circosConfig = createCircosConfig(cmd, commonConfig);
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
    private static CircosConfig createCircosConfig(CommandLine cmd, CommonConfig config) {
        return ImmutableCircosConfig.builder()
                .plotDirectory(config.outputDirectory() + File.separator + "plot")
                .circosDirectory(config.outputDirectory() + File.separator + "circos")
                .circosBinary(cmd.hasOption(CIRCOS) ? Optional.of(cmd.getOptionValue(CIRCOS)) : Optional.empty())
                .build();
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

        final String amberBaf = TumorBAFFile.generateAmberFilename(config.amberDirectory(), config.tumorSample());
        final File amberFile = new File(amberBaf);
        if (amberFile.exists()) {
            return builder.bafFile(Optional.of(amberFile)).build();
        }

        final String purpleBaf = TumorBAFFile.generatePurpleFilename(config.outputDirectory(), config.tumorSample());
        final File purpleFile = new File(purpleBaf);
        if (purpleFile.exists()) {
            return builder.bafFile(Optional.of(purpleFile)).build();
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
        throw new ParseException("Cached baf file " + purpleBaf + " not found. Please supply one of -baf or -baf_vcf arguments.");
    }

    private static void printHelp(Options opt) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Purity Ploidy Estimator (PURPLE)", opt);
    }

}
