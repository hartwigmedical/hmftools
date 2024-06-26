package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_ZIP_EXTENSION;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.drivers.DriverPrep;
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

            case DRIVER:
                return new DriverPrep(mConfig);

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

    public void extractSingleSample(boolean keepDataItems)
    {
        List<DataItem> dataItems = new ArrayList<>();

        for(CategoryType categoryType : mConfig.Categories)
        {
            CategoryPrep categoryPrep = createCategoryPrep(categoryType);
            SampleOneCategoryTask sampleTask = new SampleOneCategoryTask(0, mConfig, categoryPrep, null);
            sampleTask.run();
            dataItems.addAll(sampleTask.mDataItems);
        }

        if(keepDataItems)
            mDataItems = dataItems;

        String outputPath = getOutputPath(null);
        DataItemsIO.writeDataItemList(dataItems, outputPath);
    }

    public DataItemMatrix extractMultiSampleOneCategory(CategoryType categoryType)
    {
        CUP_LOGGER.info("Extracting category({})", categoryType);

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

    public void extractMultiSample(boolean keepDataItems)
    {
        CUP_LOGGER.info("Extracting CUPPA features in multi sample mode: {} samples, {} threads",
                mConfig.SampleIds.size(), mConfig.Threads);

        int i = 0;
        for(CategoryType categoryType : mConfig.Categories)
        {
            DataItemMatrix dataItemMatrix = extractMultiSampleOneCategory(categoryType);

            if(keepDataItems)
            {
                mDataItemMatricesByCategory.put(categoryType, dataItemMatrix);
            }

            if(mConfig.WriteByCategory)
            {
                String outputPath = getOutputPath(categoryType);
                DataItemsIO.writeDataItemMatrix(dataItemMatrix, outputPath, false);
            }
            else
            {
                String outputPath = getOutputPath(null);
                boolean append = (i != 0);
                DataItemsIO.writeDataItemMatrix(dataItemMatrix, outputPath, append);
            }
            i++;
        }
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
            extractSingleSample(keepDataItems);
        } else {
            extractMultiSample(keepDataItems);
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
