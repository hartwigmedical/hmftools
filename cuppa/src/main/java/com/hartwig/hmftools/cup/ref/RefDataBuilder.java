package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.classifierEnabled;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.feature.RefFeatures;
import com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions;
import com.hartwig.hmftools.cup.rna.RefGeneExpression;
import com.hartwig.hmftools.cup.traits.RefSampleTraits;
import com.hartwig.hmftools.cup.somatics.RefSomatics;
import com.hartwig.hmftools.cup.svs.RefSvData;

import org.jetbrains.annotations.NotNull;

public class RefDataBuilder
{
    private final RefDataConfig mConfig;

    private final SampleDataCache mSampleDataCache;

    private final List<RefClassifier> mClassifiers;

    public RefDataBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new RefDataConfig(configBuilder);

        mSampleDataCache = new SampleDataCache();

        loadSampleData(configBuilder);

        mClassifiers = Lists.newArrayList();

        // build / load traits first since some subsequent classifiers use its data (eg purity & ploidy)
        if(RefSampleTraits.requiresBuild(mConfig))
            mClassifiers.add(new RefSampleTraits(mConfig, mSampleDataCache, configBuilder));

        if(RefSomatics.requiresBuild(mConfig))
            mClassifiers.add(new RefSomatics(mConfig, mSampleDataCache, configBuilder));

        if(RefSvData.requiresBuild(mConfig))
            mClassifiers.add(new RefSvData(mConfig, mSampleDataCache));

        if(RefFeatures.requiresBuild(mConfig))
            mClassifiers.add(new RefFeatures(mConfig, mSampleDataCache, configBuilder));

        if(RefGeneExpression.requiresBuild(mConfig))
            mClassifiers.add(new RefGeneExpression(mConfig, mSampleDataCache, configBuilder));

        if(RefAltSpliceJunctions.requiresBuild(mConfig))
            mClassifiers.add(new RefAltSpliceJunctions(mConfig, mSampleDataCache));
    }

    private void loadSampleData(final ConfigBuilder configBuilder)
    {
        mSampleDataCache.loadReferenceSampleData(configBuilder.getValue(REF_SAMPLE_DATA_FILE));

        CUP_LOGGER.info("loaded {} reference samples, {} cancer types",
                mSampleDataCache.RefSampleDataList.size(), mSampleDataCache.RefCancerSampleData.size());
    }

    public void run()
    {
        if(mSampleDataCache.RefSampleDataList.isEmpty() || mSampleDataCache.RefCancerSampleData.isEmpty())
        {
            CUP_LOGGER.info("failed to load ref sample data");
            return;
        }

        CUP_LOGGER.info("CUP building ref data sets");

        for(RefClassifier classifier : mClassifiers)
        {
            if(!classifierEnabled(classifier.categoryType(), mConfig.Categories))
                continue;

            if(!classifier.buildRefDataSets())
            {
                CUP_LOGGER.error("classifier({}) ref build failed", classifier.categoryType());
                System.exit(1);
            }
        }

        mConfig.NoiseAdjustments.writeNoiseAdjustments();

        CUP_LOGGER.info("CUP ref data building complete");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        RefDataConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        RefDataBuilder refDataBuilder = new RefDataBuilder(configBuilder);
        refDataBuilder.run();
    }
}
