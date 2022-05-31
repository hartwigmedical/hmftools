package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.classifierEnabled;

import java.util.List;

import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.feature.RefFeatures;
import com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions;
import com.hartwig.hmftools.cup.rna.RefGeneExpression;
import com.hartwig.hmftools.cup.traits.RefSampleTraits;
import com.hartwig.hmftools.cup.somatics.RefSomatics;
import com.hartwig.hmftools.cup.svs.RefSvData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class RefDataBuilder
{
    private final RefDataConfig mConfig;

    private final SampleDataCache mSampleDataCache;

    private final List<RefClassifier> mClassifiers;

    public RefDataBuilder(final CommandLine cmd)
    {
        mConfig = new RefDataConfig(cmd);

        mSampleDataCache = new SampleDataCache();

        loadSampleData(cmd);


        mClassifiers = Lists.newArrayList();

        // build / load traits first since some subsequent classifiers use its data (eg purity & ploidy)
        if(RefSampleTraits.requiresBuild(mConfig))
            mClassifiers.add(new RefSampleTraits(mConfig, mSampleDataCache, cmd));

        if(RefSomatics.requiresBuild(mConfig))
            mClassifiers.add(new RefSomatics(mConfig, mSampleDataCache, cmd));

        if(RefSvData.requiresBuild(mConfig))
            mClassifiers.add(new RefSvData(mConfig, mSampleDataCache));

        if(RefFeatures.requiresBuild(mConfig))
            mClassifiers.add(new RefFeatures(mConfig, mSampleDataCache, cmd));

        if(RefGeneExpression.requiresBuild(mConfig))
            mClassifiers.add(new RefGeneExpression(mConfig, mSampleDataCache, cmd));

        if(RefAltSpliceJunctions.requiresBuild(mConfig))
            mClassifiers.add(new RefAltSpliceJunctions(mConfig, mSampleDataCache, cmd));
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadReferenceSampleData(cmd.getOptionValue(REF_SAMPLE_DATA_FILE));

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

            classifier.buildRefDataSets();
        }

        mConfig.NoiseAdjustments.writeNoiseAdjustments();

        CUP_LOGGER.info("CUP ref data building complete");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        RefDataConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        RefDataBuilder refDataBuilder = new RefDataBuilder(cmd);
        refDataBuilder.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
