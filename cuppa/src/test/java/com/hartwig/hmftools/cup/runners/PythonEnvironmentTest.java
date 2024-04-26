package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import static org.junit.Assert.assertTrue;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Ignore;
import org.junit.Test;

public class PythonEnvironmentTest
{
    private static final PythonEnvironment ENV = new PythonEnvironment(
            PythonEnvironment.DEFAULT_PYTHON_VERSION,
            "/Users/lnguyen/Desktop/pycuppa_venv_test/"
    );

    private static final String PYCUPPA_DIR = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/src/main/python/pycuppa";

    public PythonEnvironmentTest()
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
    }

    @Ignore
    @Test
    public void canInitializeEnvironment()
    {
        boolean overwrite = true;
        if(overwrite)
        {
            ENV.removeExistingPyenv();
            ENV.removeExistingVirtualEnv();
        }

        ENV.installPyenv(true);
        ENV.installPython(true);
        ENV.createVirtualEnvironment(true);

        assertTrue(ENV.pythonPath().exists());
    }

    @Ignore
    @Test
    public void canInstallPycuppa()
    {
        ENV.pipUpgrade();
        ENV.pipInstall(PYCUPPA_DIR);
    }

    @Ignore
    @Test
    public void canInitializeEnvironmentFromJar()
    {
        /*
        ## Run these commands manually
        mvn clean install -pl cuppa
        java -cp /Users/lnguyen/Hartwig/hartwigmedical/hmftools/cuppa/target/cuppa-2.1.1-jar-with-dependencies.jar com.hartwig.hmftools.cup.runners.PycuppaInstaller /Users/lnguyen/Desktop/pycuppa_venv_test
        */
        //PycuppaInstaller.main(new String[]{ "/Users/lnguyen/Desktop/pycuppa_venv_test"});
    }
}
