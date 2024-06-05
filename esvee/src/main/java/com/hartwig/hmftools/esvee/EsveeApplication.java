package com.hartwig.hmftools.esvee;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class EsveeApplication
{
    private final AssemblyConfig mConfig;

    public EsveeApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new AssemblyConfig(configBuilder);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        SV_LOGGER.info("writing to output directory({}){}",
                mConfig.OutputDir, mConfig.OutputId != null ? format(" outputId(%s)", mConfig.OutputId) : "");

        AssemblyApplication junctionProcessor = new AssemblyApplication(mConfig);

        if(!junctionProcessor.loadJunctionFiles())
        {
            SV_LOGGER.error("failed to load junction files");
            System.exit(1);
        }

        junctionProcessor.run();
        junctionProcessor.close();

        SV_LOGGER.info("Esvee complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        AssemblyConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        EsveeApplication esvee = new EsveeApplication(configBuilder);
        esvee.run();
    }
}
