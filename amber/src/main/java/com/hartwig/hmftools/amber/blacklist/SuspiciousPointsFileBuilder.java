package com.hartwig.hmftools.amber.blacklist;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SuspiciousPointsFileBuilder
{
    private final AmberTrainingConfig mConfig;

    public SuspiciousPointsFileBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new AmberTrainingConfig(configBuilder);

        if(mConfig.SampleIds.isEmpty())
        {
            AMB_LOGGER.error("no sample IDs loaded");
            System.exit(1);
        }
    }

    public void run()
    {
        AMB_LOGGER.info("running Amber suspicious points file generation from {} samples", mConfig.SampleIds.size());
        File outputFile = new File(mConfig.OutputFile);
        AmberPanelTrainer trainer = new AmberPanelTrainer(mConfig.SampleIds, new File(mConfig.AmberDir), outputFile);
        try
        {
            trainer.run();
        }
        catch(IOException e)
        {
            AMB_LOGGER.warn("Training failed", e);
        }
        AMB_LOGGER.info("Amber suspicious points file generation complete");
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        AmberTrainingConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        SuspiciousPointsFileBuilder fileBuilder = new SuspiciousPointsFileBuilder(configBuilder);
        fileBuilder.run();
    }
}
