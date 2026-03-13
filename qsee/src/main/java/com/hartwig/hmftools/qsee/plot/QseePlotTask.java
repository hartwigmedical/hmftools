package com.hartwig.hmftools.qsee.plot;

import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.RExecutor;

import org.jetbrains.annotations.Nullable;

public class QseePlotTask implements Runnable
{
    private static final String SCRIPT_PATH = "plot_qc.R";
    private static final String NO_ARG = "NA";

    private final String mVisDataPath;
    private final String mPlotPath;

    private final String mTumorId;
    @Nullable private final String mReferenceId;

    @Nullable private final String mCohortPercentilesFile;
    private final boolean mShowPlotWarnings;
    private final boolean mIsSinglePatient;

    public QseePlotTask(
            String tumorId,
            @Nullable String referenceId,
            String visDataPath,
            @Nullable String cohortPercentilesFile,
            String plotPath,
            boolean showPlotWarnings,
            boolean isSinglePatient
    ){
        mVisDataPath = visDataPath;
        mPlotPath = plotPath;

        mTumorId = tumorId;
        mReferenceId = referenceId;

        mCohortPercentilesFile = cohortPercentilesFile;
        mShowPlotWarnings = showPlotWarnings;
        mIsSinglePatient = isSinglePatient;
    }

    @Override
    public void run()
    {
        try
        {
            String referenceIdArg = (mReferenceId == null) ? NO_ARG : mReferenceId;
            String cohortPercentilesFile = (mCohortPercentilesFile == null) ? NO_ARG : mCohortPercentilesFile;
            String showPlotWarnings = Boolean.toString(mShowPlotWarnings);
            String logLevel = QC_LOGGER.getLevel().toString();
            String isSinglePatient = Boolean.toString(mIsSinglePatient);

            String[] scriptArgs = {
                    mTumorId,
                    referenceIdArg,
                    mVisDataPath,
                    cohortPercentilesFile,
                    mPlotPath,
                    showPlotWarnings,
                    isSinglePatient,
                    logLevel
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
}
