package com.hartwig.hmftools.qsee.cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QSEE_FILE_ID;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CategoryPrepFactory;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.SamplePrepTask;

public class CohortPercentilesTrainer
{
    private final TrainConfig mTrainConfig;
    private final CommonPrepConfig mCommonPrepConfig;

    private final double[] mPercentiles;

    private final List<FeaturePercentiles> mFeaturePercentiles = new ArrayList<>();

    private static final DecimalFormat PERCENTILE_FORMAT = new DecimalFormat("0.########");
    private static final DecimalFormat REF_VALUE_FORMAT = new DecimalFormat("0.########");

    public CohortPercentilesTrainer(final TrainConfig trainConfig)
    {
        mTrainConfig = trainConfig;
        mCommonPrepConfig = trainConfig.CommonPrep;

        mPercentiles = mTrainConfig.hasPercentileIncrement()
                ? PercentileTransformer.withIncrement(mTrainConfig.PercentileIncrement).getPercentiles()
                : PercentileTransformer.withNumPercentiles(mTrainConfig.NumPercentiles).getPercentiles();
    }

    private FeatureMatrix extractMultiSampleData(CategoryPrep categoryPrep, List<String> sampleIds, SampleType sampleType, String logPrefix)
    {
        FeatureMatrix sampleFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), sampleIds);

        List<Runnable> samplePrepTasks = new ArrayList<>();
        for(int sampleIndex = 0; sampleIndex < sampleIds.size(); ++sampleIndex)
        {
            SamplePrepTask task = new SamplePrepTask(
                    categoryPrep, sampleIds, sampleIndex, sampleType, sampleFeatureMatrix,
                    mCommonPrepConfig.AllowMissingInput, logPrefix
            );

            samplePrepTasks.add(task);
        }

        TaskExecutor.executeRunnables(samplePrepTasks, mCommonPrepConfig.Threads);
        samplePrepTasks.clear();

        return sampleFeatureMatrix;
    }

    private List<FeaturePercentiles> calcPercentiles(FeatureMatrix sampleFeatureMatrix, SampleType sampleType)
    {
        FeatureMatrix percentileFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), getPercentileNames());
        List<Runnable> featureTransformTasks = new ArrayList<>();

        for(int featureIndex = 0; featureIndex < sampleFeatureMatrix.numFeatures(); ++featureIndex)
        {
            PercentileTransformer transformer = PercentileTransformer.withNumPercentiles(mTrainConfig.NumPercentiles);

            double[] featureValues = sampleFeatureMatrix.getColumnValues(featureIndex);
            FeatureKey featureKey = sampleFeatureMatrix.getFeatureKeys().get(featureIndex);

            transformer.fit(featureValues, featureKey);
            percentileFeatureMatrix.addColumn(featureKey, transformer.getRefValues(), featureKey.sourceTool());
        }

        TaskExecutor.executeRunnables(featureTransformTasks, mCommonPrepConfig.Threads);
        featureTransformTasks.clear();

        Comparator<FeatureKey> comparator = Comparator.comparing(FeatureKey::type, Comparator.nullsLast(Comparator.naturalOrder()));
        percentileFeatureMatrix.getFeatureKeys().sort(comparator);

        List<FeaturePercentiles> cohortPercentiles = new ArrayList<>();
        for(FeatureKey key : percentileFeatureMatrix.getFeatureKeys())
        {
            FeaturePercentiles featurePercentiles = new FeaturePercentiles(
                    sampleType, key, mPercentiles, percentileFeatureMatrix.getColumnValues(key));

            cohortPercentiles.add(featurePercentiles);
        }

        return cohortPercentiles;
    }

    private void calcPercentilesFor(SampleType sampleType)
    {
        List<CategoryPrep> categoryPreps = new CategoryPrepFactory(mCommonPrepConfig).createCategoryPreps();
        List<String> sampleIds = mCommonPrepConfig.getSampleIds(sampleType);

        for(CategoryPrep categoryPrep : categoryPreps)
        {
            String logPrefix = String.format("sampleType(%s) category(%s) -", sampleType, categoryPrep.name());

            QC_LOGGER.info("{} Extracting sample data", logPrefix);
            FeatureMatrix sampleFeatureMatrix = extractMultiSampleData(categoryPrep, sampleIds, sampleType, logPrefix);

            QC_LOGGER.info("{} Calculating percentiles", logPrefix);
            List<FeaturePercentiles> categoryPercentiles = calcPercentiles(sampleFeatureMatrix, sampleType);

            mFeaturePercentiles.addAll(categoryPercentiles);
        }
    }

    private void writePercentileRefValues()
    {
        try(BufferedWriter writer = createBufferedWriter(getOutputFilename()))
        {
            QC_LOGGER.info("Writing percentile ref values to {}", getOutputFilename());

            StringJoiner header = new StringJoiner(TSV_DELIM);

            header.add("SampleType");
            header.add("FeatureType");
            header.add("FeatureName");
            header.add("SourceTool");
            getPercentileNames().forEach(percentileName -> header.add("Pct_" + percentileName));

            writer.write(header.toString());
            writer.newLine();

            for(FeaturePercentiles featurePercentiles : mFeaturePercentiles)
            {
                FeatureKey featureKey = featurePercentiles.key();

                FeatureType featureType = featureKey.type();
                SourceTool sourceTool = featureKey.sourceTool();
                SampleType sampleType = featurePercentiles.sampleType();

                StringJoiner line = new StringJoiner(TSV_DELIM);
                line.add(sampleType.name());
                line.add(featureType.name());
                line.add(featureKey.name());
                line.add(sourceTool.toString());

                double[] refValues = featurePercentiles.refValues();

                String refValuesStr = Arrays.stream(refValues)
                        .mapToObj(REF_VALUE_FORMAT::format)
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
        QC_LOGGER.info("Using percentiles: {}", getPercentilesString());

        if(!mCommonPrepConfig.TumorIds.isEmpty())
            calcPercentilesFor(SampleType.TUMOR);

        if(!mCommonPrepConfig.ReferenceIds.isEmpty())
            calcPercentilesFor(SampleType.REFERENCE);

        writePercentileRefValues();
    }

    private String getOutputFilename()
    {
        return mCommonPrepConfig.OutputDir + File.separator + "cohort." + QSEE_FILE_ID + ".percentiles.tsv.gz";
    }

    private List<String> getPercentileNames()
    {
        return Arrays.stream(mPercentiles).mapToObj(PERCENTILE_FORMAT::format).toList();
    }

    private String getPercentilesString()
    {
        int SHOW_ALL_THRESHOLD = 10;

        if(mPercentiles.length <= SHOW_ALL_THRESHOLD)
        {
            return Arrays.stream(mPercentiles)
                    .mapToObj(PERCENTILE_FORMAT::format)
                    .collect(Collectors.joining(", "));
        }
        else
        {
            return String.format(
                    "%s, %s, %s, ..., %s, %s, %s",
                    PERCENTILE_FORMAT.format(mPercentiles[0]),
                    PERCENTILE_FORMAT.format(mPercentiles[1]),
                    PERCENTILE_FORMAT.format(mPercentiles[2]),
                    PERCENTILE_FORMAT.format(mPercentiles[mPercentiles.length - 3]),
                    PERCENTILE_FORMAT.format(mPercentiles[mPercentiles.length - 2]),
                    PERCENTILE_FORMAT.format(mPercentiles[mPercentiles.length - 1])
            );
        }
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CommonPrepConfig.registerConfig(configBuilder);
        TrainConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        TrainConfig config = new TrainConfig(configBuilder);

        new CohortPercentilesTrainer(config).run();
    }
}

