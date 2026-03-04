package com.hartwig.hmftools.qsee.cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

public class CohortFeaturesWriter
{
    // Feature matrix is written transposed: rows as features, columns as samples

    private final QseePrepConfig mConfig;
    private final SampleType mSampleType;
    private String mOutputFile;
    private BufferedWriter mWriter;

    public CohortFeaturesWriter(QseePrepConfig config, SampleType sampleType)
    {
        mConfig = config;
        mSampleType = sampleType;

        initialise();
    }

    public void initialise()
    {
        List<String> sampleIds = mConfig.getSampleIds(mSampleType);

        mOutputFile = QseeFileCommon.generateCohortFilename(
                mConfig.OutputDir,
                "features." + mSampleType.toString().toLowerCase(),
                mConfig.OutputId,
                "tsv.gz"
        );

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile);
            StringJoiner header = new StringJoiner(TSV_DELIM);

            header.add(COL_SOURCE_TOOL);
            header.add(COL_FEATURE_TYPE);
            header.add(COL_FEATURE_NAME);
            sampleIds.forEach(header::add);

            writer.write(header.toString());
            writer.newLine();

            mWriter = writer;
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to initialise output file: {}", mOutputFile, e);
            System.exit(1);
        }
    }

    public void writeCategory(FeatureMatrix featureMatrix)
    {
        if(mWriter == null)
            throw new IllegalStateException("Feature matrix writer not yet initialised");

        try
        {
            for(int featureIndex = 0; featureIndex < featureMatrix.numFeatures(); featureIndex++)
            {
                StringJoiner line = new StringJoiner(TSV_DELIM);

                FeatureKey featureKey = featureMatrix.getFeatureKeys().get(featureIndex);
                line.add(featureKey.sourceTool().toString());
                line.add(featureKey.type().toString());
                line.add(featureKey.name());

                double[] featureValuesPerSample = featureMatrix.getColumnValues(featureIndex);
                for(double featureValue : featureValuesPerSample)
                {
                    String featureValueStr = QseeFileCommon.DECIMAL_FORMAT.format(featureValue);
                    line.add(featureValueStr);
                }

                mWriter.write(line.toString());
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to write feature matrix values to file: {}", mOutputFile, e);
            System.exit(1);
        }
    }

    public void close()
    {
        try
        {
            mWriter.close();
        }
        catch (IOException e)
        {
            QC_LOGGER.error("Failed to close feature matrix writer: {}", mOutputFile, e);
            System.exit(1);
        }
    }
}
