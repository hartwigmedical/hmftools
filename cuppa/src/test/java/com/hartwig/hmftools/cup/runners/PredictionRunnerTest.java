package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Ignore;
import org.junit.Test;

public class PredictionRunnerTest
{
    private static final String SAMPLE_ID = "TUMOR_SAMPLE";
    private static final String SAMPLE_DATA_DIR = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/src/test/resources/pipeline_output/TUMOR_SAMPLE";
    private static final String CLASSIFIER_PATH = "/Users/lnguyen/Hartwig/cloud_source_repos/common-resources-public/cuppa/37/cuppa_classifier.37.pickle.gz";
    private static final String OUTPUT_DIR = "/Users/lnguyen/Desktop/pycuppa_output/";
    private static final String VIRTUAL_ENV_PATH = "/Users/lnguyen/Desktop/pycuppa_venv_test";

    @Ignore
    @Test
    public void canPredictFromInputFeatures()
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);

        String[] args = new String[] {
                "-sample","COLO829v003T",
                "-classifier_path", CLASSIFIER_PATH,
                "-features_path", "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/src/main/python/pycuppa/resources/mock_data/input_data/new_format/COLO829v003T.cuppa_data.tsv.gz",
                "-output_dir", OUTPUT_DIR,
                "-virtual_env_path", VIRTUAL_ENV_PATH
        };

        PredictionRunner.main(args);
    }

    @Ignore
    @Test
    public void canPredictFromPipelineOutput()
    {
        String[] args = new String[] {
                "-sample",SAMPLE_ID,
                "-sample_data_dir", SAMPLE_DATA_DIR,
                "-classifier_path", CLASSIFIER_PATH,
                "-output_dir", OUTPUT_DIR,
                "-virtual_env_path", VIRTUAL_ENV_PATH
        };

        PredictionRunner.main(args);
    }
}
