package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordConstants.APP_NAME;
import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import com.hartwig.hmftools.chord.predict.ChordModel;
import com.hartwig.hmftools.chord.prep.ChordDataPrep;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class ChordApplication
{
    ChordConfig mConfig;

    public ChordApplication(ChordConfig config)
    {
        mConfig = config;
    }

    public ChordApplication(ConfigBuilder config)
    {
        mConfig = new ChordConfig(config);
    }

    public void run()
    {
        try
        {
            ChordDataPrep prep = new ChordDataPrep(mConfig);
            prep.run();

            ChordModel model;
            if(mConfig.ChordModelFile == null)
            {
                model = new ChordModel(mConfig.OutputDir + "/CHORD.tmp.rds", true);
            } else {
                model = new ChordModel(mConfig.ChordModelFile, false);
            }

            String mutContextsPath = prep.mOutputFile;
            String predictionsPath = ChordOutput.predictionsFile(mConfig);
            model.predict(mutContextsPath, predictionsPath);
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("Failed to run CHORD: " + e);
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        ChordConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ChordApplication runner = new ChordApplication(configBuilder);
        runner.run();
    }
}