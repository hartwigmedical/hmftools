package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import java.io.File;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.prep.CuppaDataPrep;
import com.hartwig.hmftools.cup.prep.PrepConfig;

import org.apache.logging.log4j.Level;

public class PredictionRunner
{
    public final PrepConfig mPrepConfig;
    public final PycuppaPredictConfig mPycuppaConfig;

    private String mFeaturesPath;

    public PredictionRunner(final ConfigBuilder configBuilder)
    {
        mPrepConfig = new PrepConfig(configBuilder);
        mPycuppaConfig = new PycuppaPredictConfig(configBuilder);
    }

    public void createOutputDirIfNotExist()
    {
        File outputDir = new File(mPrepConfig.OutputDir);
        if(!outputDir.exists())
        {
            CUP_LOGGER.info("Creating output dir: " + outputDir);
            outputDir.mkdir();
        }
    }

    public void extractFeatures()
    {
        if(!mPycuppaConfig.FeaturesPath.isEmpty())
        {
            CUP_LOGGER.info("Using pre-extracted features at: " + mPycuppaConfig.FeaturesPath);
            mFeaturesPath =  mPycuppaConfig.FeaturesPath;
            return;
        }

        CUP_LOGGER.info("Extracting features to: " + mPycuppaConfig.FeaturesPath);
        CuppaDataPrep prep = new CuppaDataPrep(mPrepConfig);
        prep.run();

        mFeaturesPath = prep.getOutputPath(null);
    }

    public void predict()
    {
        String[] command = new String[] {
                "python3 -m cuppa.predict",

                "--sample_id", mPycuppaConfig.SampleId,
                "--classifier_path", mPycuppaConfig.ClassifierPath,
                "--output_dir", mPycuppaConfig.OutputDir,

                "--features_path", mFeaturesPath
        };

        PythonEnv pythonEnvironment = new PythonEnv(PycuppaInstaller.PYTHON_VERSION, PycuppaInstaller.PYCUPPA_VENV_NAME);
        pythonEnvironment.checkRequiredPackage(PycuppaInstaller.PYCUPPA_PKG_NAME);
        new PythonEnvCommand(pythonEnvironment, String.join(" ", command)).logLevel(Level.INFO).showCommand().run();
    }

    public static void main(String[] args)
    {
        ConfigBuilder config = new ConfigBuilder(APP_NAME);

        config.disableWarnOnRepeatedRegos();
        PrepConfig.registerConfig(config);
        PycuppaPredictConfig.registerConfig(config);

        config.checkAndParseCommandLine(args);

        PredictionRunner runner = new PredictionRunner(config);
        runner.createOutputDirIfNotExist();
        runner.extractFeatures();
        runner.predict();
    }
}
