package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_ZIP_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_CATEGORY;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_KEY;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_SOURCE;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_VALUE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.feature.FeaturePrep;
import com.hartwig.hmftools.cup.rna.AltSpliceJunctionPrep;
import com.hartwig.hmftools.cup.rna.GeneExpressionPrep;
import com.hartwig.hmftools.cup.somatics.SomaticVariantPrep;
import com.hartwig.hmftools.cup.svs.StructuralVariantPrep;
import com.hartwig.hmftools.cup.traits.SampleTraitPrep;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CuppaDataPrep
{
    public final PrepConfig mConfig;

    private static final String DELIMITER = TSV_DELIM;

    @Nullable public List<DataItem> mDataItems; // only used for tests
    @Nullable public HashMap<CategoryType, DataItemMatrix> mDataItemMatricesByCategory = new HashMap<>(); // only used for tests

    public CuppaDataPrep(final ConfigBuilder configBuilder)
    {
        mConfig = new PrepConfig(configBuilder);
    }

    public CuppaDataPrep(final PrepConfig prepConfig)
    {
        mConfig = prepConfig;
    }

    public CategoryPrep createCategoryPrep(CategoryType categoryType)
    {
        switch(categoryType)
        {
            case SNV:
                return new SomaticVariantPrep(mConfig);

            case SV:
                return new StructuralVariantPrep(mConfig);

            case SAMPLE_TRAIT:
                return new SampleTraitPrep(mConfig);

            case FEATURE:
                return new FeaturePrep(mConfig);

            case ALT_SJ:
                return new AltSpliceJunctionPrep(mConfig);

            case GENE_EXP:
                return new GeneExpressionPrep(mConfig);

            default:
                throw new IllegalArgumentException(String.format("Invalid %s value", categoryType.getClass().getSimpleName()));
        }
    }

    public String getOutputPath(@Nullable CategoryType categoryType)
    {
        String path = mConfig.OutputDir + "/";

        if(mConfig.isSingleSample())
        {
            path += mConfig.SampleIds.get(0) + ".cuppa_data" + TSV_ZIP_EXTENSION;
            return path;
        }

        path += "cuppa_data.cohort";

        if(categoryType != null)
            path += "." + categoryType.toString().toLowerCase();

        if(mConfig.OutputId != null)
            path += "." + mConfig.OutputId;

        path += TSV_ZIP_EXTENSION;

        return path;
    }

    private class SingleSampleTask
    {
        public List<DataItem> extractData()
        {
            List<DataItem> dataItems = new ArrayList<>();

            for(CategoryType categoryType : mConfig.Categories)
            {
                CategoryPrep categoryPrep = createCategoryPrep(categoryType);
                SampleOneCategoryTask sampleTask = new SampleOneCategoryTask(0, mConfig, categoryPrep, null);
                sampleTask.run();
                dataItems.addAll(sampleTask.getDataItems());
            }

            return dataItems;
        }

        public void writeData(List<DataItem> dataItems, String path)
        {
            try
            {
                CUP_LOGGER.info("Writing data to: " + path);

                BufferedWriter writer = createBufferedWriter(path, false);
                StringJoiner joiner = new StringJoiner(DELIMITER);

                joiner.add(FLD_SOURCE).add(FLD_CATEGORY).add(FLD_KEY).add(FLD_VALUE);
                String header = joiner.toString();
                writer.write(header);
                writer.newLine();

                for(DataItem dataItem : dataItems)
                {
                    joiner = new StringJoiner(DELIMITER);
                    joiner.add(dataItem.Index.Source.toString());
                    joiner.add(dataItem.Index.Type.getAlias());
                    joiner.add(dataItem.Index.Key);
                    joiner.add(dataItem.Value);

                    writer.write(joiner.toString());
                    writer.newLine();
                }

                closeBufferedWriter(writer);
            }
            catch(IOException e)
            {
                CUP_LOGGER.error("Failed to write features");
                System.exit(1);
            }
        }

        public void run(boolean keepDataItems)
        {
            CUP_LOGGER.info("Extracting CUPPA features in single sample mode for sample({})", mConfig.SampleIds.get(0));

            List<DataItem> dataItems = extractData();
            if(keepDataItems)
            {
                mDataItems = dataItems;
            }

            String outputPath = getOutputPath(null);
            writeData(dataItems, outputPath);
        }

        public void run(){ run(false); }
    }

    private class MultiSampleTask
    {
        public DataItemMatrix extractDataOneCategory(CategoryType categoryType)
        {
            CUP_LOGGER.info("Extracting categoryType({})", categoryType);

            ConcurrentHashMap<DataItem.Index, String[]> featureBySampleMatrix = new ConcurrentHashMap<>();

            List<SampleOneCategoryTask> sampleTasks = new ArrayList<>();
            for(int sampleIndex = 0; sampleIndex < mConfig.SampleIds.size(); ++sampleIndex)
            {
                CategoryPrep categoryPrep = createCategoryPrep(categoryType);
                sampleTasks.add(new SampleOneCategoryTask(sampleIndex, mConfig, categoryPrep, featureBySampleMatrix));
            }

            List<Callable> callableTasks = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableTasks, mConfig.Threads);

            DataItemMatrix matrix = new DataItemMatrix(mConfig.SampleIds, featureBySampleMatrix);
            matrix.sortIndexes();

            return matrix;
        }

        public void writeDataOneCategory(DataItemMatrix dataItemMatrix, String path, boolean append)
        {
            try
            {
                CUP_LOGGER.info("Writing data to: " + path);

                StringJoiner joiner = new StringJoiner(DELIMITER);
                BufferedWriter writer = createBufferedWriter(path, append);

                if(!append)
                {
                    joiner.add(FLD_SOURCE).add(FLD_CATEGORY).add(FLD_KEY);

                    for(String sampleId : dataItemMatrix.SampleIds)
                        joiner.add(sampleId);

                    String header = joiner.toString();
                    writer.write(header);
                    writer.newLine();
                }

                for(DataItem.Index index : dataItemMatrix.getIndexes())
                {
                    joiner = new StringJoiner(DELIMITER);

                    joiner.add(index.Source.toString());
                    joiner.add(index.Type.getAlias());
                    joiner.add(index.Key);

                    for(int sampleIndex = 0; sampleIndex < dataItemMatrix.nSamples(); sampleIndex++)
                    {
                        String value = dataItemMatrix.get(index)[sampleIndex];
                        joiner.add(value);
                    }

                    writer.write(joiner.toString());
                    writer.newLine();
                }

                closeBufferedWriter(writer);
            }
            catch(Exception e)
            {
                CUP_LOGGER.error("Failed to write multi-sample feature matrix");
                System.exit(1);
            }
        }

        public void run(boolean keepDataItems)
        {
            CUP_LOGGER.info("Extracting CUPPA features in multi sample mode: {} samples, {} threads",
                    mConfig.SampleIds.size(), mConfig.Threads);

            int i = 0;
            for(CategoryType categoryType : mConfig.Categories)
            {
                DataItemMatrix dataItemMatrix = extractDataOneCategory(categoryType);

                if(keepDataItems)
                {
                    mDataItemMatricesByCategory.put(categoryType, dataItemMatrix);
                }

                if(mConfig.WriteByCategory)
                {
                    String outputPath = getOutputPath(categoryType);
                    writeDataOneCategory(dataItemMatrix, outputPath, false);
                }
                else
                {
                    String outputPath = getOutputPath(null);
                    boolean append = (i != 0);
                    writeDataOneCategory(dataItemMatrix, outputPath, append);
                }
                i++;
            }
        }

        public void run(){ run(false); }
    }

    public void run(boolean keepDataItems)
    {
        if(mConfig.SampleIds.isEmpty())
        {
            CUP_LOGGER.error("No sample ID(s) loaded");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        if(mConfig.isSingleSample())
        {
            new SingleSampleTask().run(keepDataItems);
        } else {
            new MultiSampleTask().run(keepDataItems);
        }

        CUP_LOGGER.info("Cuppa data extraction complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public void run()
    {
        run(false);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(configBuilder);
        cuppaDataPrep.run();
    }
}
