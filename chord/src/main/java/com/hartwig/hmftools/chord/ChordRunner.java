package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordConstants.APP_NAME;
import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.r.RExecutor;

public class ChordRunner
{
    ChordConfig mConfig;

    private static final String SCRIPT_RESOURCE_PATH = "./extractSigPredictHRD.R";

    public ChordRunner(ChordConfig config)
    {
        mConfig = config;
    }

    public ChordRunner(ConfigBuilder config)
    {
        mConfig = new ChordConfig(config);
    }

    public void run()
    {
        try
        {
            String sampleId = mConfig.SampleIds.get(0);

            int result = RExecutor.executeFromClasspath(
                    SCRIPT_RESOURCE_PATH,
                    true,
                    mConfig.OutputDir,
                    sampleId,
                    mConfig.purpleSomaticVcfFile(sampleId),
                    mConfig.purpleSvVcfFile(sampleId),
                    mConfig.RefGenVersion.toString()
            );

            if(result != 0)
            {
                throw new IOException("R execution failed for script: " + SCRIPT_RESOURCE_PATH);
            }
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("Failed to run CHORD: " + e);
            System.exit(1);
        }
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        ChordConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ChordRunner runner = new ChordRunner(configBuilder);
        runner.run();
    }
}