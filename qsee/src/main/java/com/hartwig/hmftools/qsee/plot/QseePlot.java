package com.hartwig.hmftools.qsee.plot;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.prep.VisDataFile;

import org.jetbrains.annotations.Nullable;

public class QseePlot
{
    private final List<String> mTumorIds;
    private final List<String> mReferenceIds;
    private final String mVisDataFile;
    @Nullable private final String mCohortPercentilesFile;
    private final String mOutputDir;
    @Nullable private final String mOutputId;
    private final boolean mShowPlotWarnings;

    private static final String SCRIPT_PATH = "plot_qc.R";

    private static final String NO_ARG = "NA";

    public QseePlot(QseePlotConfig config)
    {
        mTumorIds = config.TumorIds;
        mReferenceIds = config.ReferenceIds;
        mVisDataFile = config.VisDataFile;
        mCohortPercentilesFile = config.CohortPercentilesFile;
        mOutputDir = config.OutputDir;
        mOutputId = config.OutputId;
        mShowPlotWarnings = config.ShowPlotWarnings;
    }

    private String generateFilename(String tumorId, @Nullable String outputId)
    {
        return QseeFileCommon.generateFilename(mOutputDir, tumorId, "vis.report", outputId, "pdf");
    }

    public void plotOneSample(String tumorId, @Nullable String referenceId, @Nullable String outputId)
    {
        try
        {
            String visDataFile;
            if(mVisDataFile == null)
            {
                visDataFile = VisDataFile.generateFilename(mOutputDir, tumorId, outputId);
            }
            else
            {
                QC_LOGGER.debug("Using existing vis data file: {}", mVisDataFile);
                visDataFile = mVisDataFile;
            }

            String[] scriptArgs = {
                    tumorId,
                    referenceId == null ? NO_ARG : referenceId,
                    visDataFile,
                    mCohortPercentilesFile == null ? NO_ARG : mCohortPercentilesFile,
                    generateFilename(tumorId, outputId),
                    Boolean.toString(mShowPlotWarnings),
                    QC_LOGGER.getLevel().toString()
            };

            int exitCode = RExecutor.executeFromClasspath(SCRIPT_PATH, true, scriptArgs);

            if(exitCode != 0)
            {
                throw new RuntimeException();
            }
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

        boolean isSingleSample = mTumorIds.size() == 1;

        if(isSingleSample)
        {
            String tumorId = mTumorIds.get(0);
            String referenceId = mReferenceIds.isEmpty() ? NO_ARG : mReferenceIds.get(0);
            plotOneSample(tumorId, referenceId, mOutputId);
        }
        else
        {
            for(int sampleIndex = 0; sampleIndex < mTumorIds.size(); sampleIndex++)
            {
                String tumorId = mTumorIds.get(sampleIndex);
                String referenceId = mReferenceIds.get(sampleIndex);

                QC_LOGGER.debug("Plotting sample {}", tumorId);

                plotOneSample(tumorId, referenceId, mOutputId);
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