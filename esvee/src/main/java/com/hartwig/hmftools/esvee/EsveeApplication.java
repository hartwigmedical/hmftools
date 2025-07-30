package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.assembly.AssemblyApplication;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.caller.CallerApplication;
import com.hartwig.hmftools.esvee.caller.CallerConfig;
import com.hartwig.hmftools.esvee.depth.DepthAnnotator;
import com.hartwig.hmftools.esvee.depth.DepthConfig;
import com.hartwig.hmftools.esvee.prep.PrepApplication;
import com.hartwig.hmftools.esvee.prep.PrepConfig;

public class EsveeApplication
{
    private final ConfigBuilder mConfigBuilder;

    public EsveeApplication(final ConfigBuilder configBuilder)
    {
        mConfigBuilder = configBuilder;
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        // run prep
        runPrep();

        // run assembly
        runAssembly();

        // run depth annotation
        runDepthAnnotation();

        // run calling
        runCaller();

        SV_LOGGER.info("Esvee complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void runPrep()
    {
        PrepApplication prepApplication = new PrepApplication(mConfigBuilder);
        prepApplication.run();
        System.gc();
    }

    private void runAssembly()
    {
        AssemblyApplication assemblyApplication = new AssemblyApplication(mConfigBuilder, true);
        assemblyApplication.run();
        assemblyApplication.close();
        System.gc();
    }

    private void runDepthAnnotation()
    {
        DepthAnnotator depthAnnotator = new DepthAnnotator(mConfigBuilder);
        depthAnnotator.run();
        System.gc();
    }

    private void runCaller()
    {
        CallerApplication callerApplication = new CallerApplication(mConfigBuilder);
        callerApplication.run();
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        PrepConfig.registerConfig(configBuilder);
        AssemblyConfig.registerConfig(configBuilder);
        DepthConfig.registerConfig(configBuilder);
        CallerConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        EsveeApplication esvee = new EsveeApplication(configBuilder);
        esvee.run();
    }
}
