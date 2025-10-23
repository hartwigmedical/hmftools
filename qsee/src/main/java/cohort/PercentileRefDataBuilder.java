package cohort;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

import com.hartwig.hmftools.common.perf.TaskExecutor;

import prep.CategoryPrep;
import prep.CategoryPrepFactory;
import prep.PrepConfig;
import prep.SamplePrepTask;

public class PercentileRefDataBuilder
{
    private final PrepConfig mConfig;

    private static final int NUM_PERCENTILES = 11;
    private static final DecimalFormat PERCENTILE_FORMAT = new DecimalFormat("0.##");

    public PercentileRefDataBuilder(final PrepConfig config)
    {
        mConfig = config;
    }

    public FeatureMatrix extractSampleData()
    {
        FeatureMatrix sampleFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), mConfig.SampleIds.size());

        List<Runnable> samplePrepTasks = new ArrayList<>();
        for(int sampleIndex = 0; sampleIndex < mConfig.SampleIds.size(); ++sampleIndex)
        {
            List<CategoryPrep> categoryPreps = new CategoryPrepFactory(mConfig).createCategoryPreps();
            SamplePrepTask task = new SamplePrepTask(mConfig, sampleIndex, categoryPreps, sampleFeatureMatrix);
            samplePrepTasks.add(task);
        }

        TaskExecutor.executeRunnables(samplePrepTasks, mConfig.Threads);
        samplePrepTasks.clear();

        sampleFeatureMatrix = sampleFeatureMatrix.reorderRows(mConfig.SampleIds);

        return sampleFeatureMatrix;
    }

    public FeatureMatrix calcPercentiles(FeatureMatrix sampleFeatureMatrix)
    {
        FeatureMatrix percentileFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), NUM_PERCENTILES);

        double[] percentiles = PercentileTransformer.withNumPercentiles(NUM_PERCENTILES).getPercentiles();

        String[] percentileNames = Arrays.stream(percentiles)
                .mapToObj(x -> PERCENTILE_FORMAT.format(x))
                .toArray(String[]::new);

        percentileFeatureMatrix.setRowIds(percentileNames);

        List<Runnable> featureTransformTasks = new ArrayList<>();
        for(int featureIndex = 0; featureIndex < sampleFeatureMatrix.numFeatures(); ++featureIndex)
        {
            PercentileTransformer transformer = PercentileTransformer.withNumPercentiles(NUM_PERCENTILES);

            PercentileTransformTask task = new PercentileTransformTask(
                    featureIndex, sampleFeatureMatrix, transformer, percentileFeatureMatrix);

            featureTransformTasks.add(task);
        }

        TaskExecutor.executeRunnables(featureTransformTasks, mConfig.Threads);
        featureTransformTasks.clear();

        percentileFeatureMatrix.reorderRows(Arrays.asList(percentileNames));

        return percentileFeatureMatrix;
    }

    public void run()
    {
        FeatureMatrix sampleFeatureMatrix = extractSampleData();
        double[][] sampleFeatureValues = sampleFeatureMatrix.getValues();

        FeatureMatrix percentileFeatureMatrix = calcPercentiles(sampleFeatureMatrix);
        double[][] percentileFeatureValues = percentileFeatureMatrix.getValuesTransposed();

        boolean test = true;
    }
}
