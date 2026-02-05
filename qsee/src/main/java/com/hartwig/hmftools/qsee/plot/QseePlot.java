package com.hartwig.hmftools.qsee.plot;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.QSEE_FILE_ID;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile;
import com.hartwig.hmftools.qsee.prep.QseePrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

import org.jetbrains.annotations.Nullable;

public class QseePlot
{
    private final List<String> mTumorIds;
    private final List<String> mReferenceIds;
    private final String mSampleFeaturesFile;
    private final String mCohortPercentilesFile;
    private final String mOutputDir;

    private static final String SCRIPT_PATH = "plot_qc.R";

    private static final String MISSING_SAMPLE_ID = "NA";

    private QseePlot(List<String> tumorIds, List<String> referenceIds,
            String sampleFeaturesFile, String cohortPercentilesFile, String outputDir)
    {
        mTumorIds = tumorIds;
        mReferenceIds = referenceIds;
        mSampleFeaturesFile = sampleFeaturesFile;
        mCohortPercentilesFile = cohortPercentilesFile;
        mOutputDir = outputDir;
    }

    public QseePlot(QseePlotConfig config)
    {
        this(config.TumorIds, config.ReferenceIds, config.SampleFeaturesFile, config.CohortPercentilesFile, config.OutputDir);
    }

    public QseePlot(QseePrepConfig config)
    {
        this(
                config.CommonPrep.TumorIds,
                config.CommonPrep.ReferenceIds,
                QseePrep.formOutputFilename(config.CommonPrep),
                CohortPercentilesFile.generateFilename(config.CommonPrep.OutputDir),
                config.CommonPrep.OutputDir
        );
    }

    private String formOutputFilename(String tumorId)
    {
        return checkAddDirSeparator(mOutputDir) + tumorId + "." + QSEE_FILE_ID + ".vis.report.pdf";
    }

    public void plotOneSample(String tumorId, @Nullable String referenceId)
    {
        try
        {
            String[] scriptArgs = {
                    tumorId,
                    referenceId == null ? MISSING_SAMPLE_ID : referenceId,
                    mSampleFeaturesFile,
                    mCohortPercentilesFile,
                    formOutputFilename(tumorId),
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
            String referenceId = mReferenceIds.isEmpty() ? MISSING_SAMPLE_ID : mReferenceIds.get(0);
            plotOneSample(tumorId, referenceId);
        }
        else
        {
            for(int sampleIndex = 0; sampleIndex < mTumorIds.size(); sampleIndex++)
            {
                String tumorId = mTumorIds.get(sampleIndex);
                String referenceId = mReferenceIds.get(sampleIndex);

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