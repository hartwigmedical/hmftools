package com.hartwig.hmftools.cup.cli;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Ignore;
import org.junit.Test;

/*
Ignore these tests as building cuppa with maven might fail on other machines due to failure to install python, pyenv, pycuppa, and required
python packages
 */
@Ignore
public class PythonEnvironmentTest
{
    private static final PythonEnv ENV = new PythonEnv(
            PythonEnv.DEFAULT_PYTHON_VERSION,
            "pycuppa_venv"
    );

    private static final String PYCUPPA_DIR = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/src/main/python/pycuppa";

    public PythonEnvironmentTest()
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
    }

    @Test
    public void canInitializeEnvironment()
    {
        PythonEnv.removeExistingPyenv();
        ENV.initialize();
        assertTrue(ENV.pythonPath().exists());
    }

    @Test
    public void canInstallPycuppa()
    {
        ENV.pipInstall(PYCUPPA_DIR, true);
    }

    @Test
    public void canInitializeEnvironmentFromJar()
    {
        /*
        ## Run these commands manually
        mvn clean install -pl cuppa
        java -cp /Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/target/cuppa-2.1.1-jar-with-dependencies.jar com.hartwig.hmftools.cup.runners.PycuppaInstaller
        */
    }
}
