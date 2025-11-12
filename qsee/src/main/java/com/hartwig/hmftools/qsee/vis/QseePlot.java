package com.hartwig.hmftools.qsee.vis;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class QseePlot
{
    private final QseePlotConfig mConfig;

    private static final String SCRIPT_PATH = "plot_qc.R";

    public QseePlot(QseePlotConfig config)
    {
        mConfig = config;
    }

    public void run()
    {
        QC_LOGGER.info("Plotting QC metrics");

        try
        {
            String[] scriptArgs = {
                    mConfig.TumorId,
                    mConfig.ReferenceId,
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