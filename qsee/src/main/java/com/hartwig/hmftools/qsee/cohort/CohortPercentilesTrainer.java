package com.hartwig.hmftools.qsee.cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile.COHORT_PERCENTILES_FILE_SUFFIX;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_NAME;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_FEATURE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SOURCE_TOOL;
import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CategoryPrepFactory;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.FeaturePrep;

public class CohortPercentilesTrainer
{
    private final TrainConfig mTrainConfig;
    private final CommonPrepConfig mCommonPrepConfig;

    private final double[] mPercentiles;

    public CohortPercentilesTrainer(final TrainConfig trainConfig)
    {
        mTrainConfig = trainConfig;
        mCommonPrepConfig = trainConfig.CommonPrep;

        mPercentiles = createTransformer().getPercentiles();
    }

    private PercentileTransformer createTransformer()
    {
        return mTrainConfig.hasPercentileIncrement()
                ? PercentileTransformer.withIncrement(mTrainConfig.PercentileIncrement)
                : PercentileTransformer.withNumPercentiles(mTrainConfig.NumPercentiles);
    }

    private List<FeaturePercentiles> calcPercentiles(FeatureMatrix sampleFeatureMatrix, SampleType sampleType)
    {
        FeatureMatrix percentileFeatureMatrix = new FeatureMatrix(new HashMap<>(), getPercentileNames());

        for(int featureIndex = 0; featureIndex < sampleFeatureMatrix.numFeatures(); ++featureIndex)
        {
            PercentileTransformer transformer = createTransformer();

            double[] featureValues = sampleFeatureMatrix.getColumnValues(featureIndex);
            FeatureKey featureKey = sampleFeatureMatrix.getFeatureKeys().get(featureIndex);

            transformer.fit(featureValues, featureKey);
            percentileFeatureMatrix.addColumn(featureKey, transformer.getRefValues(), featureKey.sourceTool());
        }

        percentileFeatureMatrix.sortFeatureKeys();

        List<FeaturePercentiles> cohortPercentiles = new ArrayList<>();
        for(FeatureKey key : percentileFeatureMatrix.getFeatureKeys())
        {
            FeaturePercentiles featurePercentiles = new FeaturePercentiles(
                    sampleType, key, mPercentiles, percentileFeatureMatrix.getColumnValues(key));

            cohortPercentiles.add(featurePercentiles);
        }

        return cohortPercentiles;
    }

    private List<FeaturePercentiles> calcPercentilesFor(SampleType sampleType)
    {
        List<String> sampleIds = mCommonPrepConfig.getSampleIds(sampleType);
        if(sampleIds.isEmpty())
        {
            QC_LOGGER.info("Skipping sampleType({}) as no samples provided", sampleType.name());
            return new ArrayList<>();
        }

        List<CategoryPrep> categoryPreps = new CategoryPrepFactory(mCommonPrepConfig).createCategoryPreps();
        List<FeaturePercentiles> cohortPercentiles = new ArrayList<>();

        boolean writeCohortFeatures = mTrainConfig.WriteCohortFeatures;
        CohortFeaturesWriter cohortFeaturesWriter = null;
        if(writeCohortFeatures)
            cohortFeaturesWriter = new CohortFeaturesWriter(mCommonPrepConfig, sampleType);

        for(CategoryPrep categoryPrep : categoryPreps)
        {
            QC_LOGGER.info("Extracting cohort features - sampleType({}) category({})", sampleType, categoryPrep.name());
            FeaturePrep featurePrep = new FeaturePrep(mCommonPrepConfig);

            FeatureMatrix sampleFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), sampleIds);
            featurePrep.prepCohortCategory(categoryPrep, sampleType, sampleFeatureMatrix);

            sampleFeatureMatrix.sortFeatureKeys();

            if(writeCohortFeatures)
            {
                QC_LOGGER.info("Writing cohort features - sampleType({}) category({})", sampleType, categoryPrep.name());
                cohortFeaturesWriter.writeCategory(sampleFeatureMatrix);
            }

            QC_LOGGER.info("Calculating percentiles - sampleType({}) category({})", sampleType, categoryPrep.name());
            List<FeaturePercentiles> categoryPercentiles = calcPercentiles(sampleFeatureMatrix, sampleType);

            cohortPercentiles.addAll(categoryPercentiles);
        }

        if(writeCohortFeatures)
            cohortFeaturesWriter.close();

        return cohortPercentiles;
    }

    private void writeCohortPercentiles(List<FeaturePercentiles> cohortPercentiles, String outputFile)
    {
        try(BufferedWriter writer = createBufferedWriter(outputFile))
        {
            QC_LOGGER.info("Writing percentile ref values to {}", outputFile);

            StringJoiner header = new StringJoiner(TSV_DELIM);

            header.add(COL_SAMPLE_TYPE);
            header.add(COL_SOURCE_TOOL);
            header.add(COL_FEATURE_TYPE);
            header.add(COL_FEATURE_NAME);
            getPercentileNames().forEach(percentileName -> header.add(CohortPercentilesFile.COL_PERCENTILE_PREFIX + "_" + percentileName));

            writer.write(header.toString());
            writer.newLine();

            for(FeaturePercentiles featurePercentiles : cohortPercentiles)
            {
                FeatureKey featureKey = featurePercentiles.featureKey();

                FeatureType featureType = featureKey.type();
                SourceTool sourceTool = featureKey.sourceTool();
                SampleType sampleType = featurePercentiles.sampleType();

                StringJoiner line = new StringJoiner(TSV_DELIM);
                line.add(sampleType.toString());
                line.add(sourceTool.toString());
                line.add(featureType.toString());
                line.add(featureKey.name());

                double[] refValues = featurePercentiles.refValues();

                String refValuesStr = Arrays.stream(refValues)
                        .mapToObj(QseeFileCommon.DECIMAL_FORMAT::format)
                        .collect(Collectors.joining(TSV_DELIM));

                line.add(refValuesStr);
                writer.write(line.toString());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to calculate cohort percentiles", e);
            System.exit(1);
        }
    }

    public void run()
    {
        QC_LOGGER.info("Running {} with percentiles: {}", this.getClass().getSimpleName(), getPercentilesString());

        List<FeaturePercentiles> cohortPercentiles = new ArrayList<>();

        List<FeaturePercentiles> tumorPercentiles = calcPercentilesFor(SampleType.TUMOR);
        List<FeaturePercentiles> referencePercentiles = calcPercentilesFor(SampleType.NORMAL);

        cohortPercentiles.addAll(tumorPercentiles);
        cohortPercentiles.addAll(referencePercentiles);

        String outputFile = mCommonPrepConfig.OutputDir + File.separator + COHORT_PERCENTILES_FILE_SUFFIX;
        writeCohortPercentiles(cohortPercentiles, outputFile);
    }

    private List<String> getPercentileNames()
    {
        return Arrays.stream(mPercentiles).mapToObj(QseeFileCommon.DECIMAL_FORMAT::format).toList();
    }

    private String getPercentilesString()
    {
        int SHOW_ALL_THRESHOLD = 10;

        if(mPercentiles.length <= SHOW_ALL_THRESHOLD)
        {
            return Arrays.stream(mPercentiles)
                    .mapToObj(QseeFileCommon.DECIMAL_FORMAT::format)
                    .collect(Collectors.joining(", "));
        }
        else
        {
            return String.format(
                    "%s, %s, %s, ..., %s, %s, %s",
                    QseeFileCommon.DECIMAL_FORMAT.format(mPercentiles[0]),
                    QseeFileCommon.DECIMAL_FORMAT.format(mPercentiles[1]),
                    QseeFileCommon.DECIMAL_FORMAT.format(mPercentiles[2]),
                    QseeFileCommon.DECIMAL_FORMAT.format(mPercentiles[mPercentiles.length - 3]),
                    QseeFileCommon.DECIMAL_FORMAT.format(mPercentiles[mPercentiles.length - 2]),
                    QseeFileCommon.DECIMAL_FORMAT.format(mPercentiles[mPercentiles.length - 1])
            );
        }
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        TrainConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        TrainConfig config = new TrainConfig(configBuilder);

        new CohortPercentilesTrainer(config).run();
    }
}

