package com.hartwig.hmftools.qsee.plot;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.QSEE_FILE_ID;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.prep.QseePrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

import org.jetbrains.annotations.Nullable;

public class QseePlot
{
    private final List<String> mTumorIds;
    private final List<String> mReferenceIds;
    private final String mSampleFeaturesFile;
    @Nullable private final String mCohortPercentilesFile;
    private final String mOutputDir;
    @Nullable private final String mOutputId;

    private static final String SCRIPT_PATH = "plot_qc.R";

    private static final String NO_ARG = "NA";

    private QseePlot(List<String> tumorIds, List<String> referenceIds, String sampleFeaturesFile,
            @Nullable String cohortPercentilesFile, String outputDir, @Nullable String outputId)
    {
        mTumorIds = tumorIds;
        mReferenceIds = referenceIds;
        mSampleFeaturesFile = sampleFeaturesFile;
        mCohortPercentilesFile = cohortPercentilesFile;
        mOutputDir = outputDir;
        mOutputId = outputId;
    }

    public QseePlot(QseePlotConfig config)
    {
        this(config.TumorIds, config.ReferenceIds, config.SampleFeaturesFile, config.CohortPercentilesFile, config.OutputDir, config.OutputId);
    }

    public QseePlot(QseePrepConfig config)
    {
        this(
                config.CommonPrep.TumorIds,
                config.CommonPrep.ReferenceIds,
                QseePrep.formOutputFilename(config.CommonPrep),
                config.CohortPercentilesFile,
                config.CommonPrep.OutputDir,
                config.CommonPrep.OutputId
        );
    }

    private String formOutputFilename(String tumorId, @Nullable String outputId)
    {
        String filename = checkAddDirSeparator(mOutputDir) + tumorId + "." + QSEE_FILE_ID + ".vis.report";

        if(outputId != null)
            filename += "." + outputId;

        filename += ".pdf";

        return filename;
    }

    public void plotOneSample(String tumorId, @Nullable String referenceId, @Nullable String outputId)
    {
        try
        {
            String[] scriptArgs = {
                    tumorId,
                    referenceId == null ? NO_ARG : referenceId,
                    mSampleFeaturesFile,
                    mCohortPercentilesFile == null ? NO_ARG : mCohortPercentilesFile,
                    formOutputFilename(tumorId, outputId),
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