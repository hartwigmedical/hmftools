package com.hartwig.hmftools.cup.prep;

import static java.lang.String.format;

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
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.ConcurrentHashMap;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
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
    public final List<CategoryPrep> mDataPreparers;

    private static final String DELIMITER = TSV_DELIM;
    private static final String NULL_VALUE_STRING = "0";

    public CuppaDataPrep(final ConfigBuilder configBuilder)
    {
        mConfig = new PrepConfig(configBuilder);
        mDataPreparers = buildDataPreparers();

        if(mConfig.SampleIds.isEmpty())
        {
            CUP_LOGGER.error("No sample ID(s) loaded");
            System.exit(1);
        }
    }

    private List<CategoryPrep> buildDataPreparers()
    {
        List<CategoryPrep> preparers = Lists.newArrayList();

        for(CategoryType categoryType : mConfig.Categories)
        {
            switch(categoryType)
            {
                case SV:
                    preparers.add(new StructuralVariantPrep(mConfig));
                    break;

                case SNV:
                    preparers.add(new SomaticVariantPrep(mConfig));
                    break;

                case SAMPLE_TRAIT:
                    preparers.add(new SampleTraitPrep(mConfig));
                    break;

                case FEATURE:
                    preparers.add(new FeaturePrep(mConfig));
                    break;

                case ALT_SJ:
                    preparers.add(new AltSpliceJunctionPrep(mConfig));
                    break;

                case GENE_EXP:
                    preparers.add(new GeneExpressionPrep(mConfig));
                    break;
            }
        }

        return preparers;
    }

    public List<DataItem> getDataSingleSample()
    {
        List<DataItem> dataItems = new ArrayList<>();
        String sampleId = mConfig.SampleIds.get(0);

        for(CategoryPrep categoryPrep : mDataPreparers)
        {
            List<DataItem> categoryDataItems = null;
            try {
                CUP_LOGGER.info("  category({})", categoryPrep.categoryType());
                categoryDataItems = categoryPrep.extractSampleData(sampleId);
            } catch(Exception e) {
                CUP_LOGGER.error("Feature extraction failed for category({})", categoryPrep.categoryType());
                System.exit(1);
            }

            assert categoryDataItems != null;
            dataItems.addAll(categoryDataItems);
        }

        return dataItems;
    }

    public static void writeDataSingleSample(List<DataItem> dataItems, String path)
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

    public DataItemMatrix getDataOneCategoryMultiSample(CategoryPrep categoryPrep)
    {
        String[] sampleIds = mConfig.SampleIds.toArray(String[]::new);
        int nSamples = sampleIds.length;

        ConcurrentHashMap<DataItemIndex, String[]> featureBySampleMatrix = new ConcurrentHashMap<>();

        for(int sampleIndex = 0; sampleIndex < nSamples; sampleIndex++)
        {
            String sampleId = sampleIds[sampleIndex];

            if(sampleIndex % 100 == 0 | sampleIndex == 1)
            {
                CUP_LOGGER.info("  sampleId({}): {}/{}", sampleId, sampleIndex+1, nSamples);
            }

            List<DataItem> dataItems = categoryPrep.extractSampleData(sampleId);

            for(DataItem dataItem : dataItems)
            {
                DataItemIndex featureIndex = dataItem.Index;

                if(featureBySampleMatrix.get(featureIndex) == null)
                    featureBySampleMatrix.put(featureIndex, new String[nSamples]);

                featureBySampleMatrix.get(featureIndex)[sampleIndex] = dataItem.Value;
            }
        }

        return new DataItemMatrix(sampleIds, featureBySampleMatrix);
    }

    public synchronized static void writeDataMultiSample(DataItemMatrix dataItemMatrix, String path, boolean append)
    {
        try
        {
            CUP_LOGGER.info("Writing data to: " + path);

            StringJoiner joiner = new StringJoiner(DELIMITER);
            BufferedWriter writer = createBufferedWriter(path, append);

            if(!append)
            {
                // TODO delete old file if not append and file exists

                joiner.add(FLD_SOURCE).add(FLD_CATEGORY).add(FLD_KEY);

                for(String sampleId : dataItemMatrix.SampleIds)
                    joiner.add(sampleId);

                String header = joiner.toString();
                writer.write(header);
                writer.newLine();
            }

            for(DataItemIndex index : dataItemMatrix.getFeatureIndexes())
            {
                joiner = new StringJoiner(DELIMITER);

                joiner.add(index.Source.toString());
                joiner.add(index.Type.getAlias());
                joiner.add(index.Key);

                for(int sampleIndex = 0; sampleIndex < dataItemMatrix.nSamples(); sampleIndex++)
                {
                    String value = dataItemMatrix.get(index)[sampleIndex];
                    if(value == null)
                        value = NULL_VALUE_STRING;

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

    public String getOutputPath(@Nullable final CategoryPrep categoryPrep)
    {
        String path = mConfig.OutputDir + "cuppa_data";

        if(mConfig.isMultiSample())
            path += ".cohort";

        if(mConfig.OutputId != null)
            path += "." + mConfig.OutputId;

        if(categoryPrep != null)
            path += "." + categoryPrep.categoryType().toString().toLowerCase();

        path += TSV_ZIP_EXTENSION;

        return path;
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        if(mConfig.isSingleSample())
        {
            CUP_LOGGER.info("Extracting CUPPA features in single sample mode for sample({})", mConfig.SampleIds.get(0));

            List<DataItem> dataItems = getDataSingleSample();
            writeDataSingleSample(dataItems, getOutputPath(null));
        }
        else
        {
            CUP_LOGGER.info("Extracting CUPPA features in multi sample mode for {} samples", mConfig.SampleIds.size());

            for(int i = 0; i < mDataPreparers.size(); i++)
            {
                CategoryPrep categoryPrep = mDataPreparers.get(i);
                DataItemMatrix dataItemMatrix = getDataOneCategoryMultiSample(categoryPrep);

                if(mConfig.WriteByCategory){
                    String path = getOutputPath(categoryPrep);
                    writeDataMultiSample(dataItemMatrix, path, false);
                }
                else
                {
                    String path = getOutputPath(null);
                    boolean append = (i != 0);
                    writeDataMultiSample(dataItemMatrix, path, append);
                }
            }
        }

        CUP_LOGGER.info("Cuppa data extraction complete, mins({})", runTimeMinsStr(startTimeMs));
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
