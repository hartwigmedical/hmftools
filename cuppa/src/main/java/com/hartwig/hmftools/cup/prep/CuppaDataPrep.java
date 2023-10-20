package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import java.util.List;

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

public class CuppaDataPrep
{
    private final SampleDataWriter mSampleDataWriter;
    private final PrepConfig mConfig;

    public CuppaDataPrep(final ConfigBuilder configBuilder)
    {
        mConfig = new PrepConfig(configBuilder);

        mSampleDataWriter = new SampleDataWriter(mConfig);
    }

    public void run()
    {
        if(!mSampleDataWriter.isValid())
        {
            System.exit(1);
        }

        if(mConfig.SampleIds.isEmpty())
        {
            CUP_LOGGER.error("no sample ID(s) loaded");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        List<CategoryPrep> dataPreparers = buildDataPreparers();

        if(mConfig.isSingleSample())
        {
            String sampleId = mConfig.SampleIds.get(0);

            CUP_LOGGER.info("sample({}) extracting Cuppa data", sampleId);

            for(CategoryPrep categoryPrep : dataPreparers)
            {
                CUP_LOGGER.debug("extracting {} data", categoryPrep.categoryType());

                List<DataItem> dataItems = categoryPrep.extractSampleData(sampleId);

                if(dataItems == null)
                {
                    CUP_LOGGER.error("invalid category({}) data", categoryPrep.categoryType());
                    System.exit(1);
                }

                mSampleDataWriter.writeSampleData(sampleId, categoryPrep.categoryType(), dataItems);
            }

            CUP_LOGGER.info("sample data extraction complete");
        }
        else
        {
            CUP_LOGGER.info("extracting Cuppa data for {} samples", mConfig.SampleIds.size());

            // TODO

            CUP_LOGGER.info("Cuppa data extraction complete, mins({})", runTimeMinsStr(startTimeMs));
        }

        mSampleDataWriter.close();
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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CuppaDataPrep cuppaDataPrep = new CuppaDataPrep(configBuilder);
        cuppaDataPrep.run();
    }
}
