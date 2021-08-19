package com.hartwig.hmftools.svtools.simulation;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.utils.ConfigUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SvSimulator
{
    private final ShatteringSim mSimShattering;

    private static final Logger LOGGER = LogManager.getLogger(SvSimulator.class);

    public SvSimulator(final CommandLine cmd, final String outputDir)
    {
        final ShatteringConfig config = new ShatteringConfig(cmd);

        if(config.isValid())
            mSimShattering = new ShatteringSim(config, outputDir);
        else
            mSimShattering = null;
    }

    public void run()
    {
        LOGGER.info("starting simulations");

        if(mSimShattering != null)
        {
            mSimShattering.run();
        }

        LOGGER.info("simulations complete");
    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(OUTPUT_DIR, true, "Output directory");
        ConfigUtils.addLoggingOptions(options);
        ShatteringConfig.addCommandLineOptions(options);
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        String outputDir = parseOutputDir(cmd);

        SvSimulator simulator = new SvSimulator(cmd, outputDir);
        simulator.run();
    }

}

