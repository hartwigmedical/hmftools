package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;

import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.pileup.PileupFile;
import com.hartwig.hmftools.common.purple.baf.TumorBAF;
import com.hartwig.hmftools.common.purple.baf.TumorBAFFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class AmberApplication {

    private static final Logger LOGGER = LogManager.getLogger(AmberApplication.class);


    private static final String SAMPLE = "sample";
    private static final String REFERENCE_PILEUP = "reference";
    private static final String TUMOR_PILEUP = "tumor";
    private static final String OUTPUT_DIR = "output_dir";

    public static void main(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {
        new AmberApplication(args);
    }

    private AmberApplication(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {

        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (!cmd.hasOption(REFERENCE_PILEUP) || !cmd.hasOption(TUMOR_PILEUP) || !cmd.hasOption(OUTPUT_DIR) || !cmd.hasOption(SAMPLE)) {
            printUsageAndExit(options);
        }

        LOGGER.info("Loading reference file {}", cmd.getOptionValue(REFERENCE_PILEUP));
        final List<Pileup> normal = PileupFile.read(cmd.getOptionValue(REFERENCE_PILEUP));

        LOGGER.info("Loading reference file {}", cmd.getOptionValue(TUMOR_PILEUP));
        final List<Pileup> tumor = PileupFile.read(cmd.getOptionValue(TUMOR_PILEUP));

        LOGGER.info("Calculating Tumor BAFs");
        final AmberBAFFactory factory = new AmberBAFFactory(0.4, 0.65, 0.5, 2.0);
        final List<TumorBAF> bafs = factory.create(normal, tumor);

        final String filename = cmd.getOptionValue(OUTPUT_DIR) + File.separator + cmd.getOptionValue(SAMPLE) + ".amber.baf";
        TumorBAFFile.write(filename, bafs);

    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("AMBER", options);
        System.exit(1);
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        try {
            return parser.parse(options, args);
        } catch (ParseException e) {
            printUsageAndExit(options);
            throw e;
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(REFERENCE_PILEUP, true, "Reference Pileup");
        options.addOption(TUMOR_PILEUP, true, "Tumor Pileup");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(SAMPLE, true, "Name of sample");

        return options;
    }
}