package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcfloader.VCFFileLoader;
import com.hartwig.hmftools.patientreporter.reports.MutationalLoad;
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

    @NotNull
    private final String runDirectory;

    public static void main(final String... args) throws ParseException, IOException, HartwigException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        String runDir = cmd.getOptionValue(RUN_DIRECTORY);

        if (runDir == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-Reporter", options);
            System.exit(1);
        }

        new PatientReporterApplication(runDir).runPatientReporter();
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);

    }

    PatientReporterApplication(@NotNull final String runDirectory) {
        this.runDirectory = runDirectory;
    }

    void runPatientReporter() throws IOException, HartwigException {
        LOGGER.info("Running patient reporter on " + runDirectory);
        List<SomaticVariant> variants = VariantFilter.passOnly(
                VCFFileLoader.loadSomaticVCF(runDirectory, SOMATIC_EXTENSION));

        ConsensusRule rule = new ConsensusRule(SlicerFactory.giabHighConfidenceSlicer(),
                SlicerFactory.cpctSlicingRegionSlicer());
        variants = rule.apply(variants);
        LOGGER.info("Variants in consensus rule = " + variants.size());
        LOGGER.info("Mutational load = " + MutationalLoad.calculate(variants));
    }
}
