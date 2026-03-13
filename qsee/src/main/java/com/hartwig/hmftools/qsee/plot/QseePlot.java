package com.hartwig.hmftools.qsee.plot;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.MULTISAMPLE_SAMPLE_ID;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.prep.VisDataFile;

import org.apache.pdfbox.io.MemoryUsageSetting;
import org.apache.pdfbox.multipdf.PDFMergerUtility;
import org.jetbrains.annotations.Nullable;

public class QseePlot
{
    private final QseePlotConfig mConfig;

    private final List<String> mTumorIds;
    private final List<String> mReferenceIds;
    private final String mVisDataFile;
    private final String mOutputDir;
    @Nullable private final String mOutputId;

    public QseePlot(QseePlotConfig config)
    {
        mConfig = config;

        mTumorIds = config.TumorIds;
        mReferenceIds = config.ReferenceIds;
        mVisDataFile = config.VisDataFile;
        mOutputDir = config.OutputDir;
        mOutputId = config.OutputId;
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
                mConfig.CohortPercentilesFile,
                mConfig.ShowPlotWarnings,
                true
        );

        plotTask.run();
    }

    private void plotMultiplePatients()
    {
        List<Runnable> plotTasks = new ArrayList<>();
        List<File> plotPaths = new ArrayList<>();

        for(int sampleIndex = 0; sampleIndex < mTumorIds.size(); sampleIndex++)
        {
            String visDataPath = getVisDataPath(mTumorIds.get(sampleIndex));
            String plotPath = generateFilename(mOutputDir, mTumorIds.get(sampleIndex), mOutputId);
            plotPaths.add(new File(plotPath));

            QseePlotTask plotTask = new QseePlotTask(
                    visDataPath,
                    plotPath,
                    mTumorIds.get(sampleIndex),
                    getReferenceId(sampleIndex),
                    mConfig.CohortPercentilesFile,
                    mConfig.ShowPlotWarnings,
                    false
            );
            plotTasks.add(plotTask);
        }

        TaskExecutor.executeRunnables(plotTasks, mConfig.Threads);

        if(mConfig.MergePlots)
        {
            mergePlots(plotPaths);
        }
    }

    private void mergePlots(List<File> plotPaths)
    {
        String destinationPath = generateFilename(mOutputDir, MULTISAMPLE_SAMPLE_ID, mOutputId);

        try
        {
            QC_LOGGER.info("Merging plots to file: {}", destinationPath);

            PDFMergerUtility merger = new PDFMergerUtility();
            merger.setDestinationFileName(destinationPath);

            for(File plotPath : plotPaths)
                merger.addSource(plotPath);

            merger.mergeDocuments(MemoryUsageSetting.setupMainMemoryOnly());
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to merge plots to file: {}", destinationPath, e);
            System.exit(1);
        }
        finally
        {
            for(File plotPath : plotPaths)
            {
                if(!plotPath.exists())
                    continue;

                QC_LOGGER.debug("Deleting plot file: {}", plotPath);
                plotPath.delete();
            }
        }
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