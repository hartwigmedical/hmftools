package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class EsveeApplication
{
    private final SvConfig mConfig;

    public EsveeApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new SvConfig(configBuilder);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        SV_LOGGER.info("writing output to VCF({}) directory({})", filenamePart(mConfig.VcfFile), mConfig.OutputDir);

        JunctionProcessor junctionProcessor = new JunctionProcessor(mConfig);

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

        SvConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        EsveeApplication esvee = new EsveeApplication(configBuilder);
        esvee.run();
    }
}
