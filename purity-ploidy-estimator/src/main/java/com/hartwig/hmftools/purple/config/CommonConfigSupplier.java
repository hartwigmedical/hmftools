package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;

import java.io.File;
import java.util.function.Supplier;

import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.copynumber.freec.FreecFileLoader;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CommonConfigSupplier implements Supplier<CommonConfig> {

    private static final Logger LOGGER = LogManager.getLogger(CommonConfig.class);

    private static final String FORCE = "force";
    private static final String FREEC_DIRECTORY = "freec_dir";
    private static final String REF_SAMPLE = "ref_sample";
    private static final String TUMOR_SAMPLE = "tumor_sample";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String OUTPUT_DIRECTORY = "output_dir";
    private static final String OUTPUT_DIRECTORY_DEFAULT = "purple";

    public static void addOptions(Options options) {
        options.addOption(REF_SAMPLE, true, "The reference sample name. Defaults to value in metadata.");
        options.addOption(TUMOR_SAMPLE, true, "The tumor sample name. Defaults to value in metadata.");
        options.addOption(RUN_DIRECTORY, true, "The path containing the data for a single run.");
        options.addOption(OUTPUT_DIRECTORY, true, "The output path. Defaults to run_dir/purple/");
        options.addOption(FREEC_DIRECTORY, true, "The freec data path. Defaults to run_dir/copyNumber/refSample_tumorSample/freec/");
        options.addOption(FORCE, false, "Force recalculation of data. Do not use cached results");
    }

    private final CommonConfig config;

    public CommonConfigSupplier(CommandLine cmd, Options opt) throws ParseException, HartwigException {
        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        if (runDirectory == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Purity Ploidy Estimator (PURPLE)", opt);
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
        config = new CommonConfig(refSample, tumorSample, outputDirectory, runDirectory, freecDirectory, cmd.hasOption(FORCE));

        LOGGER.info("Reference Sample: {}, Tumor Sample: {}", config.refSample(), config.tumorSample());
        LOGGER.info("Run Directory: {}", config.runDirectory());
        LOGGER.info("Output Directory: {}", config.outputDirectory());
    }

    @NotNull
    private static String freecDirectory(@NotNull final CommandLine cmd, @NotNull final String runDirectory,
            @NotNull final String refSample, @NotNull final String tumorSample) {
        return cmd.hasOption(FREEC_DIRECTORY)
                ? cmd.getOptionValue(FREEC_DIRECTORY)
                : FreecFileLoader.getFreecBasePath(runDirectory, refSample, tumorSample);
    }

    @Override
    public CommonConfig get() {
        return config;
    }
}
