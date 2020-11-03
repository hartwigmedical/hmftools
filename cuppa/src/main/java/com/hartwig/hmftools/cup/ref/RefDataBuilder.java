package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.REF_SAMPLE_DATA_FILE;

import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.feature.RefFeatures;
import com.hartwig.hmftools.cup.rna.RefRnaExpression;
import com.hartwig.hmftools.cup.sample.RefSampleTraits;
import com.hartwig.hmftools.cup.sigs.RefSomatics;
import com.hartwig.hmftools.cup.svs.RefSvData;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class RefDataBuilder
{
    private final RefDataConfig mConfig;

    private final SampleDataCache mSampleDataCache;

    private final RefSampleTraits mSampleTraits;
    private final RefSomatics mSomatics;
    private final RefSvData mSvAnnotation;
    private final RefFeatures mFeatures;
    private final RefRnaExpression mRnaExpression;

    public RefDataBuilder(final CommandLine cmd)
    {
        mConfig = new RefDataConfig(cmd);

        mSampleDataCache = new SampleDataCache();

        loadSampleData(cmd);

        mSampleTraits = new RefSampleTraits(mConfig, mSampleDataCache);
        mSomatics = new RefSomatics(mConfig, mSampleDataCache);
        mSvAnnotation = new RefSvData(mConfig, mSampleDataCache);
        mFeatures = new RefFeatures(mConfig, mSampleDataCache);
        mRnaExpression = new RefRnaExpression(mConfig, mSampleDataCache, cmd);
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadReferenceSampleData(cmd.getOptionValue(REF_SAMPLE_DATA_FILE), false);

        CUP_LOGGER.info("loaded {} reference samples, {} cancer types",
                mSampleDataCache.SampleIds.size(), mSampleDataCache.RefCancerSampleData.size());
    }

    public void run()
    {
        CUP_LOGGER.info("CUP building ref data sets");

        mSampleTraits.buildRefDataSets();
        mSomatics.buildRefDataSets();
        mFeatures.buildRefDataSets();
        mSvAnnotation.buildRefDataSets();
        mRnaExpression.buildRefDataSets();

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
