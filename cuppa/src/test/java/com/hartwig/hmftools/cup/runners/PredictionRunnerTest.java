package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Ignore;
import org.junit.Test;

public class PredictionRunnerTest
{
    @Ignore
    @Test
    public void callingPredictionRunnerWithoutArgsShowsUsage()
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);

        String[] args = new String[] {
                "-sample","COLO829v003T",
                "-classifier_path", "/Users/lnguyen/Hartwig/cloud_source_repos/common-resources-public/cuppa/37/cuppa_classifier.37.pickle.gz",
                "-features_path", "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/src/main/python/pycuppa/resources/mock_data/input_data/new_format/COLO829v003T.cuppa_data.tsv.gz",
                "-output_dir", "/Users/lnguyen/Desktop/pycuppa_output",
                "-virtual_env_path", "/Users/lnguyen/Desktop/pycuppa_env_test"
        };

        PredictionRunner.main(args);
    }
}
