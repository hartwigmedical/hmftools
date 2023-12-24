package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.processor.Processor;

public class EsveeApplication
{
    private final SvConfig mConfig;
    private final Context mContext;

    public EsveeApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new SvConfig(configBuilder);

        mContext = Context.create(mConfig);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        SV_LOGGER.info("starting Esvee");

        Processor processor = new Processor(mContext);

        processor.run();

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
