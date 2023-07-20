package com.hartwig.hmftools.statcalcs.cooc;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.statcalcs.common.StatsCommon.STAT_LOGGER;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.statcalcs.common.StatsCommon;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CoOccurenceCalcs
{
    private ThreeVarCoOccurence mThreeVarCoOccurence;
    private TwoVarCoOccurence mTwoVarCoOccurence;
    private SampleCountsCoOccurence mSampleCountsCoOccurence;

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        String outputDir = parseOutputDir(configBuilder);

        CoOccurenceCalcs statsRoutines = new CoOccurenceCalcs();

        statsRoutines.loadConfig(configBuilder, outputDir);
        statsRoutines.runStatistics();
    }

    private static void registerConfig(final ConfigBuilder configBuilder)
    {
        StatsCommon.registerConfig(configBuilder);
        TwoVarCoOccurence.registerConfig(configBuilder);
        ThreeVarCoOccurence.registerConfig(configBuilder);
        SampleCountsCoOccurence.registerConfig(configBuilder);
    }

    public CoOccurenceCalcs()
    {
        mThreeVarCoOccurence = null;
        mSampleCountsCoOccurence = null;
        mTwoVarCoOccurence = null;
    }

    public boolean loadConfig(final ConfigBuilder configBuilder, final String outputDir)
    {
        boolean valid = true;

        if(SampleCountsCoOccurence.hasConfig(configBuilder))
        {
            mSampleCountsCoOccurence = new SampleCountsCoOccurence(configBuilder, outputDir);
        }

        if(ThreeVarCoOccurence.hasConfig(configBuilder))
        {
            mThreeVarCoOccurence = new ThreeVarCoOccurence(configBuilder, outputDir);
        }

        if(TwoVarCoOccurence.hasConfig(configBuilder))
        {
            mTwoVarCoOccurence = new TwoVarCoOccurence(configBuilder, outputDir);
        }

        return valid;
    }

    public void runStatistics()
    {
        if(mSampleCountsCoOccurence != null)
            mSampleCountsCoOccurence.run();

        if(mThreeVarCoOccurence != null)
            mThreeVarCoOccurence.run();

        if(mTwoVarCoOccurence != null)
            mTwoVarCoOccurence.run();

        STAT_LOGGER.info("run complete");
    }
}


