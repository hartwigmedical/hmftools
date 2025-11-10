package com.hartwig.hmftools.qsee.cohort;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
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
        List<String> percentileNames = CohortPercentilesFile.getPercentileNames(mPercentiles);
        FeatureMatrix percentileFeatureMatrix = new FeatureMatrix(new HashMap<>(), percentileNames);

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

    private List<FeaturePercentiles> runFor(SampleType sampleType)
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

    public void run()
    {
        QC_LOGGER.info("Running {} with percentiles: {}", this.getClass().getSimpleName(), getPercentilesDescription());

        List<FeaturePercentiles> cohortPercentiles = new ArrayList<>();

        List<FeaturePercentiles> tumorPercentiles = runFor(SampleType.TUMOR);
        List<FeaturePercentiles> referencePercentiles = runFor(SampleType.NORMAL);

        cohortPercentiles.addAll(tumorPercentiles);
        cohortPercentiles.addAll(referencePercentiles);

        try
        {
            String outputFile = CohortPercentilesFile.generateFilename(mCommonPrepConfig.OutputDir);
            CohortPercentilesFile.write(outputFile, mPercentiles, cohortPercentiles);
        }
        catch (IOException e)
        {
            QC_LOGGER.error("Failed to write cohort percentiles file", e);
            System.exit(1);
        }
    }

    private String getPercentilesDescription()
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

