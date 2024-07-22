package com.hartwig.hmftools.cup.cli;

import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import java.io.File;
import java.util.StringJoiner;

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
        if(mPredictionConfig.FeaturesPath != null)
        {
            CUP_LOGGER.info("Using pre-extracted features at: " + mPredictionConfig.FeaturesPath);
            mFeaturesPath = mPredictionConfig.FeaturesPath;
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

        StringJoiner args = new StringJoiner(" ");

        args.add("-m cuppa.predict");
        args.add("--classifier_path").add(mPredictionConfig.ClassifierPath);
        args.add("--features_path").add(mFeaturesPath);
        args.add("--output_dir").add(mPrepConfig.OutputDir);

        if(mPredictionConfig.SampleId != null)
            args.add("--sample_id").add(mPredictionConfig.SampleId);

        if(mPredictionConfig.ClfGroup != null)
            args.add("--clf_group").add(mPredictionConfig.ClfGroup);

        if(mPredictionConfig.CvPredictionsPath != null)
            args.add("--cv_predictions_path").add(mPredictionConfig.CvPredictionsPath);

        args.add("--log_format").add(PYTHON_LOG_FORMAT);

        ShellCommand command = pythonInterpreter.command(args.toString());
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
