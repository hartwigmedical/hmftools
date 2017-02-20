package com.hartwig.hmftools.consensusrule;

import java.io.IOException;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.consensus.ConsensusRule;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFFileWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ConsensusRuleFilterApplication {

    private static final Logger LOGGER = LogManager.getLogger(ConsensusRuleFilterApplication.class);

    private static final String HIGH_CONFIDENCE_BED_ARGS_DESC = "The full path towards the high confidence bed";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String EXTREME_CONFIDENCE_BED_ARGS_DESC = "The full path towards the extreme confidence bed";
    private static final String EXTREME_CONFIDENCE_BED = "extreme_confidence_bed";

    private static final String RUN_DIRECTORY_ARGS_DESC = "The full path towards the run directory";
    private static final String RUN_DIRECTORY = "run_dir";

    private static final String INPUT_VCF_ARGS_DESC = "The full path towards the input VCF to apply consensus rule to";
    private static final String INPUT_VCF = "in_vcf";

    private static final String OUTPUT_VCF_ARGS_DESC = "The full path where the filtered VCF will be written to.";
    private static final String OUTPUT_VCF = "out_vcf";

    public static void main(final String... args)
            throws ParseException, IOException, XMLStreamException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String extremeConfidenceBed = cmd.getOptionValue(EXTREME_CONFIDENCE_BED);
        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        final String inputVcf = cmd.getOptionValue(INPUT_VCF);
        final String outputVcf = cmd.getOptionValue(OUTPUT_VCF);

        if (highConfidenceBed == null || extremeConfidenceBed == null || (runDirectory == null && inputVcf == null)
                || outputVcf == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Consensus-Rule-Filter", options);
            System.exit(1);
        }

        final Slicer highConfidenceSlicer = SlicerFactory.fromBedFile(highConfidenceBed);
        final Slicer extremeConfidenceSlicer = SlicerFactory.fromBedFile(extremeConfidenceBed);
        final ConsensusRule consensusRule = new ConsensusRule(highConfidenceSlicer, extremeConfidenceSlicer);

        final ConsensusRuleFilterApplication application = new ConsensusRuleFilterApplication(consensusRule,
                outputVcf);

        if (inputVcf == null) {
            application.runOnRunDirectory(runDirectory);
        } else {
            application.runOnInputVcf(inputVcf);
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(HIGH_CONFIDENCE_BED, true, HIGH_CONFIDENCE_BED_ARGS_DESC);
        options.addOption(EXTREME_CONFIDENCE_BED, true, EXTREME_CONFIDENCE_BED_ARGS_DESC);
        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);
        options.addOption(INPUT_VCF, true, INPUT_VCF_ARGS_DESC);
        options.addOption(OUTPUT_VCF, true, OUTPUT_VCF_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private final ConsensusRule consensusRule;
    @NotNull
    private final String outputVcf;

    private ConsensusRuleFilterApplication(@NotNull final ConsensusRule consensusRule,
            @NotNull final String outputVcf) {
        this.consensusRule = consensusRule;
        this.outputVcf = outputVcf;
    }

    private void runOnRunDirectory(@NotNull final String runDirectory) throws IOException, HartwigException {
        LOGGER.info("Loading melted input from " + runDirectory);
        processVariants(VCFFileLoader.loadSomaticVCF(runDirectory, "_melted.vcf").variants());
    }

    private void runOnInputVcf(@NotNull final String inputVcf) throws IOException, HartwigException {
        LOGGER.info("Loading explicit vcf input from " + inputVcf);
        processVariants(VCFFileLoader.loadSomaticVCF(inputVcf).variants());
    }

    private void processVariants(@NotNull final List<SomaticVariant> variants) throws IOException {
        LOGGER.info("Processing " + variants.size() + " variants in consensus rule.");

        final List<SomaticVariant> filteredVariants = consensusRule.apply(variants);
        LOGGER.info("Filtered variants on consensus rule: " + filteredVariants.size() + " variants remaining.");

        VCFFileWriter.writeSomaticVCF(outputVcf, filteredVariants);
        LOGGER.info("Written filtered variants to " + outputVcf);
    }
}
