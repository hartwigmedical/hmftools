package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.isMissense;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFFileWriter;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientreporter.slicing.SlicerFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);
    private static final String SOMATIC_EXTENSION = "_melted.vcf";

    private static final String RUN_DIRECTORY_ARGS_DESC = "A path towards a single rundir.";
    private static final String RUN_DIRECTORY = "rundir";

    private static final String CPCT_SLICING_BED_ARGS_DESC = "A path towards the CPCT slicing bed.";
    private static final String CPCT_SLICING_BED = "cpct_slicing_bed";

    private static final String HIGH_CONFIDENCE_BED_ARGS_DESC = "A path towards the high confidence bed.";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String HMF_SLICING_BED_ARGS_DESC = "A path towards the HMF slicing bed.";
    private static final String HMF_SLICING_BED = "hmf_slicing_bed";

    private static final String VCF_OUTPUT_PATH_ARGS_DESC = "A path where, if provided, vcfs will be written to.";
    private static final String VCF_OUTPUT_PATH = "vcf_output_path";

    @NotNull
    private final String runDirectory;
    @NotNull
    private final String cpctSlicingBed;
    @NotNull
    private final String highConfidenceBed;
    @NotNull
    private final String hmfSlicingBed;
    @Nullable
    private final String vcfOutputPath;

    public static void main(final String... args) throws ParseException, IOException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);
        final String cpctSlicingBed = cmd.getOptionValue(CPCT_SLICING_BED);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String hmfSlicingBed = cmd.getOptionValue(HMF_SLICING_BED);
        final String vcfOutputPath = cmd.getOptionValue(VCF_OUTPUT_PATH);

        if (runDir == null || cpctSlicingBed == null || highConfidenceBed == null || hmfSlicingBed == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-Reporter", options);
            System.exit(1);
        }

        if (vcfOutputPath != null) {
            final Path vcfPath = new File(vcfOutputPath).toPath();
            if (!Files.exists(vcfPath) || !Files.isDirectory(vcfPath)) {
                LOGGER.warn("vcf_output_path has to be an existing directory!");
                System.exit(1);
            }
        }

        new PatientReporterApplication(runDir, cpctSlicingBed, highConfidenceBed, hmfSlicingBed, vcfOutputPath).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);
        options.addOption(CPCT_SLICING_BED, true, CPCT_SLICING_BED_ARGS_DESC);
        options.addOption(HIGH_CONFIDENCE_BED, true, HIGH_CONFIDENCE_BED_ARGS_DESC);
        options.addOption(HMF_SLICING_BED, true, HMF_SLICING_BED_ARGS_DESC);
        options.addOption(VCF_OUTPUT_PATH, true, VCF_OUTPUT_PATH_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    PatientReporterApplication(@NotNull final String runDirectory, @NotNull final String cpctSlicingBed,
            @NotNull final String highConfidenceBed, @NotNull final String hmfSlicingBed,
            @Nullable final String vcfOutputPath) {
        this.runDirectory = runDirectory;
        this.cpctSlicingBed = cpctSlicingBed;
        this.highConfidenceBed = highConfidenceBed;
        this.hmfSlicingBed = hmfSlicingBed;
        this.vcfOutputPath = vcfOutputPath;
    }

    void run() throws IOException, HartwigException {
        LOGGER.info("Running patient reporter on " + runDirectory);

        final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(runDirectory, SOMATIC_EXTENSION);
        LOGGER.info(" - Extracted variants for sample " + variantFile.sample());

        final List<SomaticVariant> allVariants = variantFile.variants();
        LOGGER.info(" - Total number of variants: " + allVariants.size());

        final List<SomaticVariant> allPassedVariants = passOnly(allVariants);
        LOGGER.info(" - Number of variants after applying pass-only filter: " + allPassedVariants.size());

        final ConsensusRule consensus = new ConsensusRule(SlicerFactory.fromBedFile(highConfidenceBed),
                SlicerFactory.fromBedFile(cpctSlicingBed));
        final List<SomaticVariant> consensusPassedVariants = consensus.apply(allPassedVariants);
        LOGGER.info(" - Number of variants after applying consensus rule = " + consensusPassedVariants.size());

        final List<SomaticVariant> consensusMissenseVariants = filter(consensusPassedVariants, isMissense());
        LOGGER.info(" - Mutational load: " + consensusMissenseVariants.size());
        if (vcfOutputPath != null) {
            final String mutationalLoadVCF = variantFile.sample() + "_mutational_load_variants.vcf";
            VCFFileWriter.writeSomaticVCF(vcfOutputPath + File.separator + mutationalLoadVCF,
                    consensusMissenseVariants);
        }

        final ConsequenceRule consequence = new ConsequenceRule(SlicerFactory.fromBedFile(hmfSlicingBed));
        final List<SomaticVariant> consequencePassedVariants = consequence.apply(consensusPassedVariants);
        LOGGER.info(" - Number of consequential variants to report: " + consequencePassedVariants.size());
        if (vcfOutputPath != null) {
            final String variantsToReport = variantFile.sample() + "_variants_to_report.vcf";
            VCFFileWriter.writeSomaticVCF(vcfOutputPath + File.separator + variantsToReport,
                    consequencePassedVariants);
        }
    }
}
