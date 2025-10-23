package cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import static common.QseeConstants.QC_LOGGER;
import static common.QseeConstants.QSEE_FILE_ID;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.ConcurrentHashMap;

import com.hartwig.hmftools.common.perf.TaskExecutor;

import feature.FeatureKey;
import prep.CategoryPrep;
import prep.CategoryPrepFactory;
import prep.PrepConfig;
import prep.SamplePrepTask;

public class PercentileRefDataBuilder
{
    private final PrepConfig mConfig;

    private static final int NUM_PERCENTILES = 11;
    private static final DecimalFormat PERCENTILE_FORMAT = new DecimalFormat("0.##");
    private static final DecimalFormat REF_VALUE_FORMAT = new DecimalFormat("0.####");

    public PercentileRefDataBuilder(final PrepConfig config)
    {
        mConfig = config;
    }

    private FeatureMatrix extractSampleData()
    {
        QC_LOGGER.info("Extracting sample data");

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

    private FeatureMatrix calcPercentiles(FeatureMatrix sampleFeatureMatrix)
    {
        QC_LOGGER.info("Transforming feature values to percentiles");

        FeatureMatrix percentileFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), NUM_PERCENTILES);

        double[] percentiles = PercentileTransformer.withNumPercentiles(NUM_PERCENTILES).getPercentiles();

        List<String> percentileNames = Arrays.stream(percentiles)
                .mapToObj(x -> PERCENTILE_FORMAT.format(x))
                .toList();

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

        percentileFeatureMatrix.reorderRows(percentileNames);

        return percentileFeatureMatrix;
    }

    private void writeToFile(FeatureMatrix percentileFeatureMatrix)
    {
        try
        {
            String outputFile = mConfig.OutputDir + File.separator + "cohort." + QSEE_FILE_ID + ".percentiles.tsv.gz";
            QC_LOGGER.info("Writing cohort percentile data to: {}", outputFile);

            BufferedWriter writer = createBufferedWriter(outputFile);

            // Write header
            StringJoiner header = new StringJoiner(TSV_DELIM);
            header.add("FeatureType").add("FeatureName");

            List<String> percentileNames = percentileFeatureMatrix.getRowIds();
            for(String percentileName : percentileNames)
            {
                String percentileColumnName = "Pct_" + percentileName;
                header.add(percentileColumnName);
            }

            writer.write(header.toString());
            writer.newLine();

            // Write contents
            for(int featureIndex = 0; featureIndex < percentileFeatureMatrix.numFeatures(); featureIndex++)
            {
                FeatureKey featureKey = percentileFeatureMatrix.getFeatureKeys().get(featureIndex);

                StringJoiner line = new StringJoiner(TSV_DELIM);

                line.add(featureKey.type().toString());
                line.add(featureKey.name());

                double[] percentileValues = percentileFeatureMatrix.getColumnValues(featureIndex);

                String percentileValuesStr = Arrays.stream(percentileValues)
                        .mapToObj(REF_VALUE_FORMAT::format)
                        .collect(java.util.stream.Collectors.joining(TSV_DELIM));

                line.add(percentileValuesStr);

                writer.write(line.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to write percentile data to file: {}", e.toString());
            e.printStackTrace();
        }
    }

    public void run()
    {
        FeatureMatrix sampleFeatureMatrix = extractSampleData();

        FeatureMatrix percentileFeatureMatrix = calcPercentiles(sampleFeatureMatrix);

        writeToFile(percentileFeatureMatrix);
    }
}
