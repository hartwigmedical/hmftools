package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.LOG_DEBUG;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.createCommandLine;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.createOptions;

import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.BatchRunData;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class BachelorApplication {

    private final GermlineVcfParser mGermlineVcfParser;
    private final VariantEnricher mVariantEnricher;

    private final BachelorConfig mConfig;

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);

    private BachelorApplication(final CommandLine cmd)
    {
        mConfig = new BachelorConfig(cmd);

        mGermlineVcfParser = new GermlineVcfParser(mConfig, cmd);
        mVariantEnricher = new VariantEnricher(mConfig, cmd);
    }

    private void run()
    {
        if(!mConfig.isValid())
        {
            LOGGER.error("exiting due to invalid config");
            return;
        }

        if(!mConfig.IsBatchMode)
        {
            mGermlineVcfParser.run(mConfig.GermlineVcf, mConfig.SampleId, mConfig.OutputDir);

            final List<BachelorGermlineVariant> bachelorRecords = mGermlineVcfParser.getBachelorRecords();
            mVariantEnricher.run(bachelorRecords);
        }
        else
        {
            // to be completed

        }

        mVariantEnricher.close();

        LOGGER.info("bachelor run complete");
    }

    public static void main(final String... args)
    {
        final Options options = createOptions();

        try
        {
            final CommandLine cmd = createCommandLine(options, args);

            if (cmd.hasOption(LOG_DEBUG))
                Configurator.setRootLevel(Level.DEBUG);

            BachelorApplication bachelorApp = new BachelorApplication(cmd);
            bachelorApp.run();
        }
        catch (Exception e)
        {
            LOGGER.error("config errors: {}", e.toString());
            e.printStackTrace();
        }
    }
}
