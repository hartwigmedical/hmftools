package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.prep.CuppaDataPrep;
import com.hartwig.hmftools.cup.prep.PrepConfig;

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
        mPycuppaExecutor.runBashCommandInVirtualEnv("python3 -m cuppa.predict", Level.ERROR);
    }

    public static void main(String[] args)
    {
        args = new String[] {
                "-sample","COLO829v003T",
                "-classifier_path", "/Users/lnguyen/Hartwig/cloud_source_repos/common-resources-public/cuppa/37/cuppa_classifier.37.pickle.gz",
                "-features_path", "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/src/main/python/pycuppa/resources/mock_data/input_data/new_format/COLO829v003T.cuppa_data.tsv.gz",
                "-output_dir", "/Users/lnguyen/Desktop/pycuppa_output",
                "-virtual_env_path", "/Users/lnguyen/Desktop/pycuppa_env_test"
        };

        ConfigBuilder config = new ConfigBuilder(APP_NAME);
        PredictionConfig.registerConfig(config);
        config.checkAndParseCommandLine(args);

        PredictionRunner runner = new PredictionRunner(config);
        runner.predict();
    }
}
