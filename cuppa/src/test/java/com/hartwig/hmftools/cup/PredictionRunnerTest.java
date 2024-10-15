package com.hartwig.hmftools.cup;

import static java.lang.String.format;

import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.cup.cli.PredictionRunner;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.After;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class PredictionRunnerTest
{
    private static final String CUPPA_DIR = System.getProperty("user.dir");
    private static final String USER_HOME = System.getProperty("user.home");

    private static final String SAMPLE_ID = "TUMOR_SAMPLE";
    private static final String SAMPLE_DATA_DIR = CUPPA_DIR + "/src/test/resources/pipeline_output/TUMOR_SAMPLE";
    private static final String CLASSIFIER_PATH = format("%s/Hartwig/common-resources-public/cuppa/37/cuppa_classifier.37.pickle.gz", USER_HOME);
    private static final String OUTPUT_DIR = System.getProperty("java.io.tmpdir") + "/cuppa_output/";
    private static final String PYTHON_PATH = format("%s/.pyenv/versions/3.9.4/envs/pycuppa_venv/bin/python", USER_HOME);

    public PredictionRunnerTest()
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
    }

    @After
    public void teardown() throws IOException {
        FileUtils.deleteDirectory(new File(OUTPUT_DIR));
    }

    @Test
    public void canPredictFromInputFeatures()
    {
        String[] args = new String[] {
                "-sample","COLO829v003T",
                "-classifier_path", CLASSIFIER_PATH,
                "-features_path", CUPPA_DIR + "/src/main/python/pycuppa/resources/mock_data/input_data/COLO829v003T.cuppa_data.tsv.gz",
                "-output_dir", OUTPUT_DIR,
                "-python_path", PYTHON_PATH
        };

        PredictionRunner.main(args);
    }

    @Test
    public void canPredictFromPipelineOutput()
    {
        String[] args = new String[] {
                "-sample",SAMPLE_ID,
                "-sample_data_dir", SAMPLE_DATA_DIR,
                "-classifier_path", CLASSIFIER_PATH,
                "-output_dir", OUTPUT_DIR,
                "-python_path", PYTHON_PATH
        };

        PredictionRunner.main(args);
    }
}

/*
java -cp /Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/target/cuppa-2.1.1-jar-with-dependencies.jar \
com.hartwig.hmftools.cup.cli.PredictionRunner \
-sample TUMOR_SAMPLE \
-output_dir /Users/lnguyen/Desktop/pycuppa_test_output \
-classifier_path /Users/lnguyen/Hartwig/cloud_source_repos/common-resources-public/cuppa/37/cuppa_classifier.37.pickle.gz \
-sample_data_dir /Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/src/test/resources/pipeline_output/TUMOR_SAMPLE \
-python_path /Users/lnguyen/.pyenv/versions/3.9.4/envs/pycuppa_venv/bin/python
*/
