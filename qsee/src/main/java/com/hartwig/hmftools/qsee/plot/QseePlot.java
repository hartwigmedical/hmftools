package com.hartwig.hmftools.qsee.plot;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public class QseePlot
{
    private final QseePlotConfig mConfig;

    private static final String SCRIPT_PATH = "plot_qc.R";

    private static final String MISSING_SAMPLE_ID = "NA";

    public QseePlot(QseePlotConfig config)
    {
        mConfig = config;
    }

    public void plotOneSample(String tumorId, @Nullable String referenceId)
    {
        try
        {
            String[] scriptArgs = {
                    tumorId,
                    referenceId == null ? MISSING_SAMPLE_ID : referenceId,
                    mConfig.SampleFeaturesFile,
                    mConfig.CohortPercentilesFile,
                    mConfig.OutputPath,
                    QC_LOGGER.getLevel().toString()
            };

            RExecutor.executeFromClasspath(SCRIPT_PATH, true, scriptArgs);
        }
        catch(IOException | InterruptedException e)
        {
            QC_LOGGER.error("Failed to run QC plot script", e);
            System.exit(1);
        }
    }

    public void run()
    {
        QC_LOGGER.info("Plotting QC metrics");

        boolean isSingleSample = mConfig.TumorIds.size() == 1;

        if(isSingleSample)
        {
            String tumorId = mConfig.TumorIds.get(0);
            String referenceId = mConfig.ReferenceIds.isEmpty() ? MISSING_SAMPLE_ID : mConfig.ReferenceIds.get(0);
            plotOneSample(tumorId, referenceId);
        }
        else
        {
            for(int sampleIndex = 0; sampleIndex < mConfig.TumorIds.size(); sampleIndex++)
            {
                String tumorId = mConfig.TumorIds.get(sampleIndex);
                String referenceId = mConfig.ReferenceIds.get(sampleIndex);

                QC_LOGGER.debug("Plotting sample {}", tumorId);

                plotOneSample(tumorId, referenceId);
            }
        }

        QC_LOGGER.info("Completed plotting QC metrics");
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        QseePlotConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        QseePlotConfig config = new QseePlotConfig(configBuilder);

        new QseePlot(config).run();
    }
}