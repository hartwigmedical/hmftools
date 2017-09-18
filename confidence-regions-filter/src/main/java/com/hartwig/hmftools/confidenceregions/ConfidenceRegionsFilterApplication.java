package com.hartwig.hmftools.confidenceregions;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.confidenceregions.ConfidenceRegionsRule;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFFileWriter;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ConfidenceRegionsFilterApplication {

    private static final Logger LOGGER = LogManager.getLogger(ConfidenceRegionsFilterApplication.class);
    private static final String SOMATIC_EXTENSION = "_post_processed.vcf";
    private static final String FLAG_HEADER =
            "##FILTER=<ID=" + ConfidenceRegionsRule.FILTER + ",Description=\"HMF-filter: Not included in confidence regions\">";

    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String EXTREME_CONFIDENCE_BED = "extreme_confidence_bed";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String INPUT_VCF = "in_vcf";
    private static final String OUTPUT_VCF = "out_vcf";
    private static final String USE_FILTER_FLAG = "use_filter_flag";
    private static final String BATCH_MODE = "batch_mode";
    private static final String BATCH_MODE_IN_DIRECTORY = "batch_mode_in_dir";
    private static final String BATCH_MODE_OUT_DIRECTORY = "batch_mode_out_dir";

    public static void main(final String... args) throws ParseException, IOException, XMLStreamException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String extremeConfidenceBed = cmd.getOptionValue(EXTREME_CONFIDENCE_BED);
        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        final String inputVcf = cmd.getOptionValue(INPUT_VCF);
        final String outputVcf = cmd.getOptionValue(OUTPUT_VCF);
        final boolean useFilterFlag = cmd.hasOption(USE_FILTER_FLAG);
        final boolean batchMode = cmd.hasOption(BATCH_MODE);
        final String batchModeInDirectory = cmd.getOptionValue(BATCH_MODE_IN_DIRECTORY);
        final String batchModeOutDirectory = cmd.getOptionValue(BATCH_MODE_OUT_DIRECTORY);

        final boolean validBatchMode = batchMode && batchModeInDirectory != null && batchModeOutDirectory != null;
        final boolean validSingleMode = !batchMode && outputVcf != null && (inputVcf != null || runDirectory != null);
        if (highConfidenceBed == null || extremeConfidenceBed == null || !(validBatchMode || validSingleMode)) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Confidence-Regions-Rule-Filter", options);
            System.exit(1);
        }

        final Slicer highConfidenceSlicer = SlicerFactory.fromBedFile(highConfidenceBed);
        final Slicer extremeConfidenceSlicer = SlicerFactory.fromBedFile(extremeConfidenceBed);
        final ConfidenceRegionsRule confidenceRegionsRule =
                ConfidenceRegionsRule.fromSlicers(highConfidenceSlicer, extremeConfidenceSlicer);

        final ConfidenceRegionsFilterApplication application = new ConfidenceRegionsFilterApplication(confidenceRegionsRule, useFilterFlag);

        if (batchMode) {
            LOGGER.info("Running confidence regions rule filter in batch mode");
            application.runBatchModeOnDirectory(batchModeInDirectory, batchModeOutDirectory);
        } else if (inputVcf == null) {
            application.runOnRunDirectory(runDirectory, outputVcf);
        } else {
            application.runOnInputVcf(inputVcf, outputVcf);
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(HIGH_CONFIDENCE_BED, true, "The full path towards the high confidence bed");
        options.addOption(EXTREME_CONFIDENCE_BED, true, "The full path towards the extreme confidence bed");
        options.addOption(RUN_DIRECTORY, true, "The full path towards the run directory");
        options.addOption(INPUT_VCF, true, "The full path towards the input VCF to apply confidence regions rule to");
        options.addOption(OUTPUT_VCF, true, "The full path where the filtered VCF will be written to");
        options.addOption(USE_FILTER_FLAG, false,
                "If set, updates the filter flag rather than removing the variants that dont pass confidence regions rule");
        options.addOption(BATCH_MODE, false, "If set, runs the confidence regions filter in batch mode on " + BATCH_MODE_IN_DIRECTORY);
        options.addOption(BATCH_MODE_IN_DIRECTORY, true, "When in batch mode, assumes this folder contains individual runs");
        options.addOption(BATCH_MODE_OUT_DIRECTORY, true, "When in batch mode, writes files to this directory");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private final ConfidenceRegionsRule confidenceRegionsRule;
    private final boolean useFilterFlag;

    private ConfidenceRegionsFilterApplication(@NotNull final ConfidenceRegionsRule confidenceRegionsRule, final boolean useFilterFlag) {
        this.confidenceRegionsRule = confidenceRegionsRule;
        this.useFilterFlag = useFilterFlag;
    }

    private void runBatchModeOnDirectory(@NotNull final String runDirectory, @NotNull final String outputDirectory)
            throws IOException, HartwigException {
        for (final Path run : Files.list(new File(runDirectory).toPath()).collect(Collectors.toList())) {
            LOGGER.info("Processing " + run.toFile().getName());
            final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(run.toFile().getPath(), SOMATIC_EXTENSION);
            final String extension = useFilterFlag ? "_confidenceregions_flagged.vcf" : "_confidenceregions_filtered.vcf";
            final String outputVcf = outputDirectory + File.separator + variantFile.sample() + extension;
            processVariants(variantFile, outputVcf);
        }
    }

    private void runOnRunDirectory(@NotNull final String runDirectory, @NotNull final String outputVcf)
            throws IOException, HartwigException {
        LOGGER.info("Loading post processed input from " + runDirectory);
        processVariants(VCFFileLoader.loadSomaticVCF(runDirectory, SOMATIC_EXTENSION), outputVcf);
    }

    private void runOnInputVcf(@NotNull final String inputVcf, @NotNull final String outputVcf) throws IOException, HartwigException {
        LOGGER.info("Loading explicit vcf input from " + inputVcf);
        processVariants(VCFFileLoader.loadSomaticVCF(inputVcf), outputVcf);
    }

    private void processVariants(@NotNull final VCFSomaticFile inputFile, @NotNull final String outputVcf) throws IOException {
        final List<SomaticVariant> variants = inputFile.variants();
        LOGGER.info("Processing " + variants.size() + " variants in confidence regions rule for " + inputFile.sample());

        final List<SomaticVariant> filteredVariants = useFilterFlag
                ? confidenceRegionsRule.updateFilterFlagForUnreliableVariants(variants)
                : confidenceRegionsRule.removeUnreliableVariants(variants);
        LOGGER.info("Filtered variants on confidence regions rule: " + VariantFilter.passOnly(filteredVariants).size()
                + " variants remaining.");

        final List<String> additionalMetaInformation = useFilterFlag ? Lists.newArrayList(FLAG_HEADER) : Lists.newArrayList();
        VCFFileWriter.writeSomaticVCF(outputVcf, inputFile, filteredVariants, additionalMetaInformation);
        LOGGER.info("Written filtered variants to " + outputVcf);
    }
}
