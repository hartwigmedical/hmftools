package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.qc.AmberQC;
import com.hartwig.hmftools.common.amber.qc.AmberQCFactory;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.pileup.PileupFile;
import com.hartwig.hmftools.common.version.VersionInfo;

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
    private static final String MIN_HET_AF_PERCENTAGE = "min_het_af_percent";
    private static final String MAX_HET_AF_PERCENTAGE = "max_het_af_percent";
    private static final String MIN_DEPTH_PERCENTAGE = "min_depth_percent";
    private static final String MAX_DEPTH_PERCENTAGE = "max_depth_percent";

    private static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    private static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;
    private static final double DEFAULT_MIN_DEPTH_PERCENTAGE = 0.5;
    private static final double DEFAULT_MAX_DEPTH_PERCENTAGE = 1.5;

    public static void main(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {
        new AmberApplication(args);
    }

    private AmberApplication(final String... args) throws ParseException, IOException, ExecutionException, InterruptedException {
        final VersionInfo versionInfo = new VersionInfo("amber.version");
        LOGGER.info("AMBER version: {}", versionInfo.version());
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);
        if (!cmd.hasOption(REFERENCE_PILEUP) || !cmd.hasOption(TUMOR_PILEUP) || !cmd.hasOption(OUTPUT_DIR) || !cmd.hasOption(SAMPLE)) {
            printUsageAndExit(options);
        }

        final File outputDir = new File(cmd.getOptionValue(OUTPUT_DIR));
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Unable to write directory " + cmd.getOptionValue(OUTPUT_DIR));
        }

        final AmberBAFFactory factory = new AmberBAFFactory(defaultValue(cmd, MIN_HET_AF_PERCENTAGE, DEFAULT_MIN_HET_AF_PERCENTAGE),
                defaultValue(cmd, MAX_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE),
                defaultValue(cmd, MIN_DEPTH_PERCENTAGE, DEFAULT_MIN_DEPTH_PERCENTAGE),
                defaultValue(cmd, MAX_DEPTH_PERCENTAGE, DEFAULT_MAX_DEPTH_PERCENTAGE));

        LOGGER.info("Loading tumor file {}", cmd.getOptionValue(TUMOR_PILEUP));
        final List<Pileup> tumor = PileupFile.read(cmd.getOptionValue(TUMOR_PILEUP));

        LOGGER.info("Loading reference file {}", cmd.getOptionValue(REFERENCE_PILEUP));
        final List<Pileup> normal = PileupFile.read(cmd.getOptionValue(REFERENCE_PILEUP));

        LOGGER.info("Calculating BAFs");
        final List<AmberBAF> result = factory.create(normal, tumor);

        LOGGER.info("Generating QC Stats");
        final AmberQC qcStats = AmberQCFactory.create(result);
        final String qcFilename = AmberQCFile.generateFilename(cmd.getOptionValue(OUTPUT_DIR), cmd.getOptionValue(SAMPLE));

        final String filename = AmberBAFFile.generateAmberFilename(cmd.getOptionValue(OUTPUT_DIR), cmd.getOptionValue(SAMPLE));
        LOGGER.info("Persisting file {}", filename);
        AmberBAFFile.write(filename, result);
        AmberQCFile.write(qcFilename, qcStats);
        versionInfo.write(outputDir.toString());

        LOGGER.info("Complete");
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

    private static double defaultValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue) {
        if (cmd.hasOption(opt)) {
            final double result = Double.valueOf(cmd.getOptionValue(opt));
            LOGGER.info("Using non default value {} for parameter {}", result, opt);
            return result;
        }

        return defaultValue;
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(REFERENCE_PILEUP, true, "Reference Pileup");
        options.addOption(TUMOR_PILEUP, true, "Tumor Pileup");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(SAMPLE, true, "Name of sample");
        options.addOption(MIN_HET_AF_PERCENTAGE, true, "Min heterozygous AF%");
        options.addOption(MAX_HET_AF_PERCENTAGE, true, "Max heterozygous AF%");
        options.addOption(MIN_DEPTH_PERCENTAGE, true, "Max percentage of median depth");
        options.addOption(MAX_DEPTH_PERCENTAGE, true, "Min percentage of median depth");
        return options;
    }
}