package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.isMissense;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFFileWriter;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientreporter.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.util.ConsequenceCount;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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

    private static final String BATCH_MODE_ARGS_DESC = "If set, runs in batch mode (Caution!!! Korneel Only)";
    private static final String BATCH_MODE = "batch_mode";

    @NotNull
    private final String runDirectory;
    @NotNull
    private final ConsensusRule consensusRule;
    @NotNull
    private final ConsequenceRule consequenceRule;
    @Nullable
    private final String vcfOutputPath;
    private final boolean batchMode;

    public static void main(final String... args) throws ParseException, IOException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);
        final String cpctSlicingBed = cmd.getOptionValue(CPCT_SLICING_BED);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String hmfSlicingBed = cmd.getOptionValue(HMF_SLICING_BED);
        final String vcfOutputPath = cmd.getOptionValue(VCF_OUTPUT_PATH);
        final boolean batchMode = cmd.hasOption(BATCH_MODE);

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

        final ConsensusRule consensusRule = new ConsensusRule(SlicerFactory.fromBedFile(highConfidenceBed),
                SlicerFactory.fromBedFile(cpctSlicingBed));
        final ConsequenceRule consequenceRule = new ConsequenceRule(SlicerFactory.fromBedFile(hmfSlicingBed));
        new PatientReporterApplication(runDir, consensusRule, consequenceRule, vcfOutputPath, batchMode).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);
        options.addOption(CPCT_SLICING_BED, true, CPCT_SLICING_BED_ARGS_DESC);
        options.addOption(HIGH_CONFIDENCE_BED, true, HIGH_CONFIDENCE_BED_ARGS_DESC);
        options.addOption(HMF_SLICING_BED, true, HMF_SLICING_BED_ARGS_DESC);
        options.addOption(VCF_OUTPUT_PATH, true, VCF_OUTPUT_PATH_ARGS_DESC);
        options.addOption(BATCH_MODE, false, BATCH_MODE_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    PatientReporterApplication(@NotNull final String runDirectory, @NotNull final ConsensusRule consensusRule,
            @NotNull final ConsequenceRule consequenceRule, @Nullable final String vcfOutputPath,
            final boolean batchMode) {
        this.runDirectory = runDirectory;
        this.consensusRule = consensusRule;
        this.consequenceRule = consequenceRule;
        this.vcfOutputPath = vcfOutputPath;
        this.batchMode = batchMode;
    }

    void run() throws IOException, HartwigException {
        if (batchMode) {
            batchRun();
        } else {
            patientRun();
        }
    }

    private void batchRun() throws IOException, HartwigException {
        // KODU: We assume "run directory" is a path with a lot of directories on which we can all run in patient mode.
        VariantConsequence[] consequences = VariantConsequence.values();

        String header = "SAMPLE,VARIANT_COUNT,PASS_ONLY_COUNT,CONSENSUS_COUNT,MISSENSE_COUNT,CONSEQUENCE_COUNT";
        for (VariantConsequence consequence : consequences) {
            header += ("," + consequence.name() + "_COUNT");
        }
        System.out.println(header);

        for (Path run : Files.list(new File(runDirectory).toPath()).collect(Collectors.toList())) {
            final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(run.toFile().getPath(), SOMATIC_EXTENSION);

            final ReportContext report = ReportContext.fromVCFFile(variantFile, consensusRule, consequenceRule);
            final Map<VariantConsequence, Integer> counts = ConsequenceCount.count(report.consensusPassedVariants);
            String consequenceList = Strings.EMPTY;
            for (final VariantConsequence consequence : consequences) {
                consequenceList += ("," + counts.get(consequence).toString());
            }

            final String out =
                    variantFile.sample() + "," + report.allVariants.size() + "," + report.allPassedVariants.size()
                            + "," + report.consensusPassedVariants.size() + "," + report.missenseVariants.size() + ","
                            + report.consequencePassedVariants.size() + consequenceList;
            System.out.println(out);
        }
    }

    private void patientRun() throws IOException, HartwigException {
        LOGGER.info("Running patient reporter on " + runDirectory);

        final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(runDirectory, SOMATIC_EXTENSION);
        LOGGER.info("  Extracted variants for sample " + variantFile.sample());

        final ReportContext report = ReportContext.fromVCFFile(variantFile, consensusRule, consequenceRule);
        LOGGER.info("  Total number of variants: " + report.allVariants.size());
        LOGGER.info("  Number of variants after applying pass-only filter: " + report.allPassedVariants.size());
        LOGGER.info("  Number of variants after applying consensus rule = " + report.consensusPassedVariants.size());
        LOGGER.info("  Number of missense variants in consensus rule (mutational load) : "
                + report.missenseVariants.size());
        LOGGER.info("  Number of consequential variants to report: " + report.consequencePassedVariants.size());

        if (vcfOutputPath != null) {
            final String missenseVCF =
                    vcfOutputPath + File.separator + variantFile.sample() + "_missense_variants.vcf";
            VCFFileWriter.writeSomaticVCF(missenseVCF, report.missenseVariants);
            LOGGER.info("    Written missense variants to " + missenseVCF);

            final String consequenceVCF =
                    vcfOutputPath + File.separator + variantFile.sample() + "_consequential_variants.vcf";
            VCFFileWriter.writeSomaticVCF(consequenceVCF, report.consequencePassedVariants);
            LOGGER.info("    Written consequential variants to " + consequenceVCF);
        }
    }

    private static class ReportContext {
        @NotNull
        private final List<SomaticVariant> allVariants;
        @NotNull
        private final List<SomaticVariant> allPassedVariants;
        @NotNull
        private final List<SomaticVariant> consensusPassedVariants;
        @NotNull
        private final List<SomaticVariant> missenseVariants;
        @NotNull
        private final List<SomaticVariant> consequencePassedVariants;

        @NotNull
        static ReportContext fromVCFFile(@NotNull final VCFSomaticFile variantFile,
                @NotNull final ConsensusRule consensusRule, @NotNull final ConsequenceRule consequenceRule) {
            final List<SomaticVariant> allVariants = variantFile.variants();
            final List<SomaticVariant> allPassedVariants = passOnly(allVariants);
            final List<SomaticVariant> consensusPassedVariants = consensusRule.apply(allPassedVariants);
            final List<SomaticVariant> missenseVariants = filter(consensusPassedVariants, isMissense());
            final List<SomaticVariant> consequencePassedVariants = consequenceRule.apply(consensusPassedVariants);
            return new ReportContext(allVariants, allPassedVariants, consensusPassedVariants, missenseVariants,
                    consequencePassedVariants);
        }

        private ReportContext(@NotNull final List<SomaticVariant> allVariants,
                @NotNull final List<SomaticVariant> allPassedVariants,
                @NotNull final List<SomaticVariant> consensusPassedVariants,
                @NotNull final List<SomaticVariant> missenseVariants,
                @NotNull final List<SomaticVariant> consequencePassedVariants) {
            this.allVariants = allVariants;
            this.allPassedVariants = allPassedVariants;
            this.consensusPassedVariants = consensusPassedVariants;
            this.missenseVariants = missenseVariants;
            this.consequencePassedVariants = consequencePassedVariants;
        }
    }
}
