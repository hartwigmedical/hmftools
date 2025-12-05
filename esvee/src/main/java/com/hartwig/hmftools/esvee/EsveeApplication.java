package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
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
    private final List<Stage> mStages;

    private enum Stage
    {
        PREP,
        ASSEMBLY,
        REF_DEPTH,
        CALLING;
    }

    private static final String CFG_STAGE = "stages";

    public EsveeApplication(final ConfigBuilder configBuilder)
    {
        mConfigBuilder = configBuilder;

        mStages = Lists.newArrayList();

        if(configBuilder.hasValue(CFG_STAGE))
        {
            Arrays.stream(configBuilder.getValue(CFG_STAGE).split(ITEM_DELIM, -1)).forEach(x -> mStages.add(Stage.valueOf(x)));
        }
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        // run prep
        runPrep();
        System.gc();

        // run assembly
        runAssembly();
        System.gc();

        // run depth annotation
        runDepthAnnotation();
        System.gc();

        // run calling
        runCaller();

        SV_LOGGER.info("Esvee complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private boolean runStage(final Stage stage) { return mStages.isEmpty() || mStages.contains(stage); }

    private void runPrep()
    {
        if(!runStage(Stage.PREP))
            return;

        PrepApplication prepApplication = new PrepApplication(mConfigBuilder);
        prepApplication.run();
    }

    private void runAssembly()
    {
        if(!runStage(Stage.ASSEMBLY))
            return;

        AssemblyApplication assemblyApplication = new AssemblyApplication(mConfigBuilder, true);
        assemblyApplication.run();
        assemblyApplication.close();
    }

    private void runDepthAnnotation()
    {
        if(!runStage(Stage.REF_DEPTH))
            return;

        DepthAnnotator depthAnnotator = new DepthAnnotator(mConfigBuilder);
        depthAnnotator.run();
    }

    private void runCaller()
    {
        if(!runStage(Stage.CALLING))
            return;

        CallerApplication callerApplication = new CallerApplication(mConfigBuilder);
        callerApplication.run();
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(CFG_STAGE, false, enumValueSelectionAsStr(Stage.values(), "Stages"));

        PrepConfig.registerConfig(configBuilder);
        AssemblyConfig.registerConfig(configBuilder);
        DepthConfig.registerConfig(configBuilder);
        CallerConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        EsveeApplication esvee = new EsveeApplication(configBuilder);
        esvee.run();
    }
}
