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

    public PercentileRefDataBuilder(final PrepConfig config)
    {
        mConfig = config;
    }

    public FeatureMatrix extractSampleData()
    {
        FeatureMatrix sampleFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), mConfig.SampleIds.size());

        List<Runnable> sampleTasks = new ArrayList<>();
        for(int sampleIndex = 0; sampleIndex < mConfig.SampleIds.size(); ++sampleIndex)
        {
            List<CategoryPrep> categoryPreps = new CategoryPrepFactory(mConfig).createCategoryPreps();
            //List<CategoryPrep> categoryPreps = List.of(new CobaltGcMediansPrep(mConfig));
            sampleTasks.add(new SamplePrepTask(mConfig, sampleIndex, categoryPreps, sampleFeatureMatrix));
        }

        TaskExecutor.executeRunnables(sampleTasks, mConfig.Threads);

        return sampleFeatureMatrix;
    }

    private static final int NUM_PERCENTILES = 11;
    private static final DecimalFormat PERCENTILE_FORMAT = new DecimalFormat("0.##");

    public FeatureMatrix calcPercentiles(FeatureMatrix sampleFeatureMatrix)
    {
        FeatureMatrix percentileFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), NUM_PERCENTILES);

        double[] percentiles = PercentileTransformer.withNumPercentiles(NUM_PERCENTILES).getPercentiles();

        String[] percentileNames = Arrays.stream(percentiles)
                .mapToObj(x -> PERCENTILE_FORMAT.format(x))
                .toArray(String[]::new);

        percentileFeatureMatrix.setRowIds(percentileNames);

        double[][] featureSampleValues = sampleFeatureMatrix.getValuesTransposed();
        for(int featureIndex = 0; featureIndex < featureSampleValues.length; ++featureIndex)
        {
            double[] featureValues = featureSampleValues[featureIndex];

            PercentileTransformer transformer = PercentileTransformer.withNumPercentiles(NUM_PERCENTILES);
            transformer.fit(featureValues);

            String featureKey = sampleFeatureMatrix.getFeatureKeys().get(featureIndex);
            percentileFeatureMatrix.addColumn(featureKey, transformer.getRefValues());
        }
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

//    public static void main(String[] args)
//    {
//        PercentileTransformer transformer = new PercentileTransformer(11);
//        boolean test = true;
//    }
}
