package cohort;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static common.QseeConstants.APP_NAME;
import static common.QseeConstants.QC_LOGGER;
import static common.QseeConstants.QSEE_FILE_ID;

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

import common.SampleType;
import feature.FeatureKey;
import feature.FeatureType;
import feature.SourceTool;
import prep.CategoryPrep;
import prep.CategoryPrepFactory;
import prep.PrepConfig;
import prep.SamplePrepTask;

public class CohortPercentilesTrainer
{
    private final PrepConfig mConfig;

    private static final int NUM_PERCENTILES = 11;
    private static final DecimalFormat PERCENTILE_FORMAT = new DecimalFormat("0.##");
    private static final DecimalFormat REF_VALUE_FORMAT = new DecimalFormat("0.########");

    public CohortPercentilesTrainer(final PrepConfig config)
    {
        mConfig = config;
    }

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return basePath + File.separator + sampleId + "." + QSEE_FILE_ID + ".percentiles.tsv.gz";
    }

    private static List<String> getPercentileNames()
    {
        double[] percentiles = PercentileTransformer.withNumPercentiles(NUM_PERCENTILES).getPercentiles();

        return Arrays.stream(percentiles)
                .mapToObj(x -> PERCENTILE_FORMAT.format(x))
                .toList();
    }

    private void writeHeader(BufferedWriter writer) throws IOException
    {
        StringJoiner header = new StringJoiner(TSV_DELIM);
        header.add("SampleType");
        header.add("FeatureType");
        header.add("FeatureName");
        header.add("SourceTool");

        for(String percentileName : getPercentileNames())
        {
            String percentileColumnName = "Pct_" + percentileName;
            header.add(percentileColumnName);
        }

        writer.write(header.toString());
        writer.newLine();
    }

    private FeatureMatrix extractMultiSampleData(CategoryPrep categoryPrep, List<String> sampleIds, SampleType sampleType)
    {
        FeatureMatrix sampleFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), sampleIds);

        List<Runnable> samplePrepTasks = new ArrayList<>();
        for(int sampleIndex = 0; sampleIndex < sampleIds.size(); ++sampleIndex)
        {
            SamplePrepTask task = new SamplePrepTask(categoryPrep, sampleIds, sampleIndex, sampleType, sampleFeatureMatrix);
            samplePrepTasks.add(task);
        }

        TaskExecutor.executeRunnables(samplePrepTasks, mConfig.Threads);
        samplePrepTasks.clear();

        return sampleFeatureMatrix;
    }

    private void calcPercentiles(FeatureMatrix sampleFeatureMatrix, FeatureMatrix percentileFeatureMatrix)
    {
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

        Comparator<FeatureKey> comparator = Comparator.comparing(FeatureKey::type, Comparator.nullsLast(Comparator.naturalOrder()));
        percentileFeatureMatrix.getFeatureKeys().sort(comparator);
    }

    private void writePercentiles(BufferedWriter writer, FeatureMatrix percentileFeatureMatrix, SampleType sampleType) throws IOException
    {
        for(int featureIndex = 0; featureIndex < percentileFeatureMatrix.numFeatures(); featureIndex++)
        {
            FeatureKey featureKey = percentileFeatureMatrix.getFeatureKeys().get(featureIndex);
            FeatureType featureType = featureKey.type();
            SourceTool sourceTool = percentileFeatureMatrix.getSourceTool(featureKey);

            StringJoiner line = new StringJoiner(TSV_DELIM);

            line.add(sampleType.name());
            line.add(featureType.name());
            line.add(featureKey.name());
            line.add(sourceTool.toString());

            double[] percentileValues = percentileFeatureMatrix.getColumnValues(featureIndex);

            String percentileValuesStr = Arrays.stream(percentileValues)
                    .mapToObj(REF_VALUE_FORMAT::format)
                    .collect(Collectors.joining(TSV_DELIM));

            line.add(percentileValuesStr);

            writer.write(line.toString());
            writer.newLine();
        }
    }

    private void runFor(SampleType sampleType, BufferedWriter writer) throws IOException
    {
        QC_LOGGER.info("Calculating {} cohort percentiles", sampleType.toString());

        List<CategoryPrep> categoryPreps = new CategoryPrepFactory(mConfig).createCategoryPreps();
        List<String> sampleIds = mConfig.getSampleIds(sampleType);

        FeatureMatrix percentileFeatureMatrix = new FeatureMatrix(new ConcurrentHashMap<>(), getPercentileNames());
        for(CategoryPrep categoryPrep : categoryPreps)
        {
            QC_LOGGER.info("Running prep for category: {}", categoryPrep.getClass().getSimpleName());
            FeatureMatrix sampleFeatureMatrix = extractMultiSampleData(categoryPrep, sampleIds, sampleType);

            QC_LOGGER.info("Calculating percentiles for category: {}", categoryPrep.getClass().getSimpleName());
            calcPercentiles(sampleFeatureMatrix, percentileFeatureMatrix);
        }

        QC_LOGGER.info("Writing cohort percentile data to file");
        writePercentiles(writer, percentileFeatureMatrix, sampleType);

        QC_LOGGER.info("Completed calculating {} cohort percentiles", sampleType.toString());
    }

    public void run()
    {
        String outputFile = generateFilename(mConfig.OutputDir, "cohort");
        QC_LOGGER.info("Writing cohort percentile data to: {}", outputFile);

        try(BufferedWriter writer = createBufferedWriter(outputFile))
        {
            writeHeader(writer);

            if(!mConfig.TumorIds.isEmpty())
            {
                runFor(SampleType.TUMOR, writer);
            }

            if(!mConfig.ReferenceIds.isEmpty())
            {
                runFor(SampleType.REFERENCE, writer);
            }
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to write percentile data", e);
        }
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PrepConfig prepConfig = new PrepConfig(configBuilder);

        new CohortPercentilesTrainer(prepConfig).run();
    }
}

