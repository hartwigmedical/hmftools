package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.Level;

public class PredictionRunner
{
    public final PredictionConfig mConfig;
    public final PycuppaExecutor mPycuppaExecutor;

    public PredictionRunner(final PredictionConfig predictionConfig)
    {
        mConfig = predictionConfig;

        mPycuppaExecutor = new PycuppaExecutor(mConfig.VirtualEnvPath);
        mPycuppaExecutor.initialize();
    }

    public PredictionRunner(final ConfigBuilder configBuilder)
    {
        mConfig = new PredictionConfig(configBuilder);

        mPycuppaExecutor = new PycuppaExecutor(mConfig.VirtualEnvPath);
        mPycuppaExecutor.initialize();
    }

    public void predict()
    {
        String[] command = new String[] {
                "python3 -m cuppa.predict",
                "--sample_id", mConfig.SampleId,
                "--classifier_path", mConfig.ClassifierPath,
                "--features_path", mConfig.FeaturesPath,
                "--output_dir", mConfig.OutputDir,
        };

        mPycuppaExecutor.runBashCommandInVirtualEnv(String.join(" ", command), Level.INFO);
    }

    public void printUsage()
    {
        // TODO: Add usage for java component
        mPycuppaExecutor.runBashCommandInVirtualEnv("python3 -m cuppa.predict", Level.OFF);
    }

    public static void main(String[] args)
    {
        ConfigBuilder config = new ConfigBuilder(APP_NAME);
        PredictionConfig.registerConfig(config);
        config.checkAndParseCommandLine(args);

        PredictionRunner runner = new PredictionRunner(config);
        runner.predict();
    }
}
