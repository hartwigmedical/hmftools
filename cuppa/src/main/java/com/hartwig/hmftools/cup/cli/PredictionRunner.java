package com.hartwig.hmftools.cup.cli;

import static com.hartwig.hmftools.cup.utils.CuppaConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import java.io.File;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.prep.CuppaDataPrep;
import com.hartwig.hmftools.cup.prep.PrepConfig;

import org.apache.logging.log4j.Level;

public class PredictionRunner
{
    public final PrepConfig mPrepConfig;
    public final PredictionConfig mPredictionConfig;

    private String mFeaturesPath;

    private static final String PYTHON_LOG_FORMAT = "'[python] %(levelname)s %(name)s | %(message)s'";
    private static final String PYCUPPA_PKG_NAME = "pycuppa";

    public PredictionRunner(final ConfigBuilder configBuilder)
    {
        mPrepConfig = new PrepConfig(configBuilder);
        mPredictionConfig = new PredictionConfig(configBuilder);
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
        if(!mPredictionConfig.FeaturesPath.isEmpty())
        {
            CUP_LOGGER.info("Using pre-extracted features at: " + mPredictionConfig.FeaturesPath);
            mFeaturesPath =  mPredictionConfig.FeaturesPath;
            return;
        }

        CuppaDataPrep prep = new CuppaDataPrep(mPrepConfig);
        prep.run();

        mFeaturesPath = prep.getOutputPath(null);
    }

    public void predict()
    {
        PythonInterpreter pythonInterpreter = new PythonInterpreter(mPredictionConfig.PythonPath)
                .requirePackages(PYCUPPA_PKG_NAME);

        ShellCommand command = pythonInterpreter.command(
                "-m cuppa.predict",
                "--sample_id", mPredictionConfig.SampleId,
                "--classifier_path", mPredictionConfig.ClassifierPath,
                "--output_dir", mPredictionConfig.OutputDir,
                "--features_path", mFeaturesPath,
                "--log_format", PYTHON_LOG_FORMAT
        );
        command.logLevel(Level.INFO);
        CUP_LOGGER.info("Predicting using command: {}", command);
        command.run();
    }

    public static void main(String[] args)
    {
        ConfigBuilder config = new ConfigBuilder(APP_NAME);

        config.disableWarnOnRepeatedRegos();
        PrepConfig.registerConfig(config);
        PredictionConfig.registerConfig(config);

        config.checkAndParseCommandLine(args);

        PredictionRunner runner = new PredictionRunner(config);
        runner.createOutputDirIfNotExist();
        runner.extractFeatures();
        runner.predict();
    }
}
