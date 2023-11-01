package com.hartwig.hmftools.svtools.simulation;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.svtools.simulation.ShatteringConfig.registerConfig;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SvSimulator
{
    private final ShatteringSim mSimShattering;

    private static final Logger LOGGER = LogManager.getLogger(SvSimulator.class);

    public SvSimulator(final ConfigBuilder configBuilder)
    {
        String outputDir = parseOutputDir(configBuilder);

        final ShatteringConfig config = new ShatteringConfig(configBuilder);

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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        registerConfig(configBuilder);
        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        SvSimulator simulator = new SvSimulator(configBuilder);
        simulator.run();
    }
}

