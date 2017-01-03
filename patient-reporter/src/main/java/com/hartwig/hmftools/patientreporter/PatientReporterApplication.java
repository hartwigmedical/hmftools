package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.isMissense;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.vcfloader.VCFFileLoader;
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

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);
    private static final String SOMATIC_EXTENSION = "_melted.vcf";

    private static final String RUN_DIRECTORY_ARGS_DESC = "A path towards a single rundir.";
    private static final String RUN_DIRECTORY = "rundir";

    private static final String CPCT_SLICING_BED_ARGS_DESC = "A path towards the CPCT slicing bed.";
    private static final String CPCT_SLICING_BED = "cpct_slicing_bed";

    private static final String HIGH_CONFIDENCE_BED_ARGS_DESC = "A path towards the high confidence bed.";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    @NotNull
    private final String runDirectory;
    @NotNull
    private final String cpctSlicingBed;
    @NotNull
    private final String highConfidenceBed;

    public static void main(final String... args) throws ParseException, IOException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);
        final String cpctSlicingBed = cmd.getOptionValue(CPCT_SLICING_BED);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);

        if (runDir == null || cpctSlicingBed == null || highConfidenceBed == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-Reporter", options);
            System.exit(1);
        }

        new PatientReporterApplication(runDir, cpctSlicingBed, highConfidenceBed).runPatientReporter();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);
        options.addOption(CPCT_SLICING_BED, true, CPCT_SLICING_BED_ARGS_DESC);
        options.addOption(HIGH_CONFIDENCE_BED, true, HIGH_CONFIDENCE_BED_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    PatientReporterApplication(@NotNull final String runDirectory, @NotNull final String cpctSlicingBed,
            @NotNull final String highConfidenceBed) {
        this.runDirectory = runDirectory;
        this.cpctSlicingBed = cpctSlicingBed;
        this.highConfidenceBed = highConfidenceBed;
    }

    void runPatientReporter() throws IOException, HartwigException {
        LOGGER.info("Running patient reporter on " + runDirectory);
        final List<SomaticVariant> allPassedVariants = passOnly(
                VCFFileLoader.loadSomaticVCF(runDirectory, SOMATIC_EXTENSION));

        final ConsensusRule rule = new ConsensusRule(SlicerFactory.fromBedFile(highConfidenceBed),
                SlicerFactory.fromBedFile(cpctSlicingBed));
        final List<SomaticVariant> consensusPassedVariants = rule.apply(allPassedVariants);
        LOGGER.info("Number of variants after applying consensus rule = " + consensusPassedVariants.size());

        final List<SomaticVariant> consensusMissenseVariants = filter(consensusPassedVariants, isMissense());
        LOGGER.info("Mutational load: " + consensusMissenseVariants.size());
    }
}
