package com.hartwig.hmftools.ctdna.utils;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.CN_SEGMENT_FILE_ID;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.SUMMARY_FILE_ID;

import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.ctdna.purity.CopyNumberProfile;
import com.hartwig.hmftools.ctdna.purity.PurityConfig;
import com.hartwig.hmftools.ctdna.purity.SampleData;

import org.apache.commons.cli.ParseException;

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

    public static void main(final String[] args) throws ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PurityConfig.addConfig(configBuilder);

        addLoggingOptions(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);

        CnPurityPlotter cnPurityPlotter = new CnPurityPlotter(configBuilder);
        cnPurityPlotter.run();

        CT_LOGGER.info("CtDNA CN fit plotting complete");
    }
}
