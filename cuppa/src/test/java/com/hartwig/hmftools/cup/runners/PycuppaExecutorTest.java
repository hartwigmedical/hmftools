package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Ignore;
import org.junit.Test;

public class PycuppaExecutorTest
{
    @Ignore
    @Test
    public void canInitializePycuppaExecutor() throws IOException
    {
        File virtualEnvPath = new File(System.getProperty("java.io.tmpdir") + "/pycuppa_venv/");

        FileUtils.deleteDirectory(virtualEnvPath);

        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
        PycuppaExecutor executor = new PycuppaExecutor(virtualEnvPath.getPath());
        executor.initialize();

        FileUtils.deleteDirectory(virtualEnvPath);
    }

    @Ignore
    @Test
    public void canInitializePycuppaExecutorFromJar()
    {
        /*
        ## Run these commands manually
        mvn clean install -pl cuppa
        java -cp /Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/target/cuppa-2.1.1-jar-with-dependencies.jar com.hartwig.hmftools.cup.pycuppa.PycuppaExecutor "/Users/lnguyen/Desktop/pycuppa_venv_test"
        */
    }
}

//String argString = new StringJoiner(" ")
//        .add("--classifier_path").add("/Users/lnguyen/Hartwig/cloud_source_repos/common-resources-public/cuppa/37/cuppa_classifier.37.pickle.gz")
//        .add("--features_path").add("/Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/src/main/python/pycuppa/resources/mock_data/input_data/new_format/COLO829v003T.cuppa_data.tsv.gz")
//        .add("--output_dir").add("/Users/lnguyen/Desktop/pycuppa_output")
//        .toString();
//
//executor.predict(argString);

