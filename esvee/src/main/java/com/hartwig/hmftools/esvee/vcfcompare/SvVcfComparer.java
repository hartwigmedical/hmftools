package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class SvVcfComparer
{
    private final CompareConfig mConfig;

    private final String mSampleId;
    private final List<CompareTask> mCompareTasks;


    public SvVcfComparer(final CompareConfig config)
    {
        mConfig = config;

        mSampleId = config.SampleId;
        mCompareTasks = config.CompareTasks;
    }

    public void run()
    {
        SV_LOGGER.info("Comparing VCFs for sample: " + mSampleId);

        if(mCompareTasks.contains(CompareTask.MATCH_BREAKENDS))
            new BreakendMatchTask(mConfig).run();

        if(mCompareTasks.contains(CompareTask.LINE_COMPARE))
            new LineCompareTask(mConfig).run();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CompareConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CompareConfig compareConfig = new CompareConfig(configBuilder);

        SvVcfComparer svVcfComparer = new SvVcfComparer(compareConfig);
        svVcfComparer.run();
    }
}
