package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.LOG_DEBUG;

import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.sample.RefSampleTraits;
import com.hartwig.hmftools.cup.sigs.RefSignatures;

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
    private final RefSignatures mSignatures;

    public RefDataBuilder(final CommandLine cmd)
    {
        mConfig = new RefDataConfig(cmd);

        mSampleDataCache = new SampleDataCache();

        loadSampleData(cmd);
        mSampleTraits = new RefSampleTraits(mConfig, mSampleDataCache);
        mSignatures = new RefSignatures(mConfig, mSampleDataCache);
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadReferenceSampleData(cmd.getOptionValue(REF_SAMPLE_DATA_FILE));
    }

    public void run()
    {
        mSampleTraits.buildRefDataSets();
        mSignatures.buildRefDataSets();

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
