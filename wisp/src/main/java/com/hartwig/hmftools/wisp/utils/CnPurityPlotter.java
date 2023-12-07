package com.hartwig.hmftools.wisp.utils;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.wisp.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.CN_SEGMENT_FILE_ID;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.SUMMARY_FILE_ID;

import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.wisp.purity.cn.CopyNumberProfile;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.SampleData;

public class CnPurityPlotter
{
    private final PurityConfig mConfig;

    public CnPurityPlotter(final ConfigBuilder configBuilder)
    {
        mConfig = new PurityConfig(configBuilder);

        if(mConfig.Samples.isEmpty())
            System.exit(1);
    }

    public void run()
    {
        String summaryFile = mConfig.formFilename(SUMMARY_FILE_ID);
        String cnSegmentsFile = mConfig.formFilename(CN_SEGMENT_FILE_ID);

        if(!Files.exists(Paths.get(summaryFile)))
        {
            CT_LOGGER.error("summary file({}) missing for plots", cnSegmentsFile);
            System.exit(1);
        }
        else if(!Files.exists(Paths.get(cnSegmentsFile)))
        {
            CT_LOGGER.error("copy number segments file({}) missing for plots", cnSegmentsFile);
            System.exit(1);
        }

        for(SampleData sample : mConfig.Samples)
        {
            for(String sampleId : sample.CtDnaSamples)
            {
                CopyNumberProfile.plotCopyNumberGcRatioFit(sample.PatientId, sampleId, mConfig);
            }
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PurityConfig.addConfig(configBuilder);

        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        CnPurityPlotter cnPurityPlotter = new CnPurityPlotter(configBuilder);
        cnPurityPlotter.run();

        CT_LOGGER.info("CtDNA CN fit plotting complete");
    }
}
