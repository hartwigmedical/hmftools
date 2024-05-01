package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class PredictionRunnerTest
{
    private static final String CUPPA_DIR = System.getProperty("user.dir");

    private static final String SAMPLE_ID = "TUMOR_SAMPLE";
    private static final String SAMPLE_DATA_DIR = CUPPA_DIR + "/src/test/resources/pipeline_output/TUMOR_SAMPLE";
    private static final String CLASSIFIER_PATH = "/Users/lnguyen/Hartwig/cloud_source_repos/common-resources-public/cuppa/37/cuppa_classifier.37.pickle.gz";
    private static final String OUTPUT_DIR = System.getProperty("java.io.tmpdir") + "/cuppa_output/";

    @Test
    public void canPredictFromInputFeatures() throws IOException
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);

        String[] args = new String[] {
                "-sample","COLO829v003T",
                "-classifier_path", CLASSIFIER_PATH,
                "-features_path", CUPPA_DIR + "/src/main/python/pycuppa/resources/mock_data/input_data/new_format/COLO829v003T.cuppa_data.tsv.gz",
                "-output_dir", OUTPUT_DIR,
        };

        PredictionRunner.main(args);

        FileUtils.deleteDirectory(new File(OUTPUT_DIR));
    }

    @Test
    public void canPredictFromPipelineOutput() throws IOException
    {
        String[] args = new String[] {
                "-sample",SAMPLE_ID,
                "-sample_data_dir", SAMPLE_DATA_DIR,
                "-classifier_path", CLASSIFIER_PATH,
                "-output_dir", OUTPUT_DIR,
        };

        PredictionRunner.main(args);

        FileUtils.deleteDirectory(new File(OUTPUT_DIR));
    }
}
