package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.LOG_DEBUG;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.createOptions;

import java.util.List;

import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

public class BachelorApplication {

    private final GermlineVcfParser mGermlineVcfParser;
    private final VariantEnricher mVariantEnricher;

    private final BachelorConfig mConfig;

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
            BACH_LOGGER.error("Exiting due to invalid config");
            return;
        }

        if(!mGermlineVcfParser.run(mConfig.GermlineVcf, mConfig.SampleId))
        {
            BACH_LOGGER.error("Germline VCF parse failed");
            return;
        }

        final List<BachelorGermlineVariant> bachelorRecords = mGermlineVcfParser.getBachelorRecords();
        mVariantEnricher.run(bachelorRecords);

        BACH_LOGGER.info("Bachelor run complete");
    }

    public static void main(final String... args)
    {
        final Options options = createOptions();

        try
        {
            final CommandLine cmd = new DefaultParser().parse(options, args);

            if (cmd.hasOption(LOG_DEBUG))
                Configurator.setRootLevel(Level.DEBUG);

            BachelorApplication bachelorApp = new BachelorApplication(cmd);
            bachelorApp.run();
        }
        catch (Exception e)
        {
            BACH_LOGGER.error("Config errors: {}", e.toString());
            e.printStackTrace();
        }
    }
}
