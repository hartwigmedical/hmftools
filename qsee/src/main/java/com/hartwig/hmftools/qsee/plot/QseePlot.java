package com.hartwig.hmftools.qsee.plot;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.MULTISAMPLE_SAMPLE_ID;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.perf.TaskExecutor;
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
    private final int mThreads ;

    public QseePlot(QseePlotConfig config)
    {
        mTumorIds = config.TumorIds;
        mReferenceIds = config.ReferenceIds;
        mVisDataFile = config.VisDataFile;
        mCohortPercentilesFile = config.CohortPercentilesFile;
        mOutputDir = config.OutputDir;
        mOutputId = config.OutputId;
        mShowPlotWarnings = config.ShowPlotWarnings;
        mThreads = config.Threads;
    }

    public static String generateFilename(String basePath, String tumorId, @Nullable String outputId)
    {
        return QseeFileCommon.generateFilename(basePath, tumorId, "vis.report", outputId, "pdf");
    }

    private boolean isSinglePatient() { return mTumorIds.size() <= 1 && mReferenceIds.size() <= 1; }

    private String getVisDataPath(String tumorId)
    {
        if(mVisDataFile != null)
        {
            QC_LOGGER.debug("Using existing vis data file: {}", mVisDataFile);
            return mVisDataFile;
        }

        return isSinglePatient()
                ? VisDataFile.generateFilename(mOutputDir, tumorId, mOutputId)
                : VisDataFile.generateFilename(mOutputDir, MULTISAMPLE_SAMPLE_ID, mOutputId);
    }

    private String getReferenceId(int sampleIndex)
    {
        boolean isTumorOnlyMode = mReferenceIds.isEmpty();
        return isTumorOnlyMode ? null : mReferenceIds.get(sampleIndex);
    }

    private void plotOnePatient()
    {
        String visDataPath = getVisDataPath(mTumorIds.get(0));
        String plotPath = generateFilename(mOutputDir, mTumorIds.get(0), mOutputId);

        QseePlotTask plotTask = new QseePlotTask(
                visDataPath,
                plotPath,
                mTumorIds.get(0),
                getReferenceId(0),
                mCohortPercentilesFile,
                mShowPlotWarnings,
                true
        );

        plotTask.run();
    }

    private void plotMultiplePatients()
    {
        List<Runnable> plotTasks = new ArrayList<>();

        for(int sampleIndex = 0; sampleIndex < mTumorIds.size(); sampleIndex++)
        {
            String visDataPath = getVisDataPath(mTumorIds.get(sampleIndex));
            String plotPath = generateFilename(mOutputDir, mTumorIds.get(sampleIndex), mOutputId);

            QseePlotTask plotTask = new QseePlotTask(
                    visDataPath,
                    plotPath,
                    mTumorIds.get(sampleIndex),
                    getReferenceId(sampleIndex),
                    mCohortPercentilesFile,
                    mShowPlotWarnings,
                    false
            );
            plotTasks.add(plotTask);
        }

        TaskExecutor.executeRunnables(plotTasks, mThreads);
    }

    public void run()
    {
        QC_LOGGER.info("Plotting QC metrics");

        if(isSinglePatient()) {
            plotOnePatient();
        } else {
            plotMultiplePatients();
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