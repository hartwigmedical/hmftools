package com.hartwig.hmftools.sage.tinc;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class TincAnalyser
{
    private final TincConfig mConfig;
    private final VariantCache mVariantCache;

    public TincAnalyser(final TincConfig config)
    {
        mConfig = config;
        mVariantCache = new VariantCache(mConfig);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        mVariantCache.loadVariants();

        TincCalculator tincCalculator = new TincCalculator(mConfig, mVariantCache.fittingVariants());
        tincCalculator.run();

        SG_LOGGER.info("TINC complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        TincConfig.registerFullConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        TincConfig config = new TincConfig(configBuilder);
        TincAnalyser tincAnalyser = new TincAnalyser(config);
        tincAnalyser.run();
    }
}
