package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.File;

import com.google.common.io.Resources;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

public class PycuppaInstaller
{
    public static final String PYTHON_VERSION = "3.9.4";

    private static final File TMP_DIR = new File(System.getProperty("java.io.tmpdir"));
    private static final File PYCUPPA_TMP_DIR = new File(TMP_DIR + "/pycuppa/");
    private static final File PYCUPPA_RESOURCE_DIR = new File(Resources.getResource("pycuppa/").getPath());

    public static final String PYCUPPA_VENV_NAME = "pycuppa_venv";
    public static final String PYCUPPA_PKG_NAME = "pycuppa";

    public final PythonEnv mPythonEnv;

    public PycuppaInstaller()
    {
        mPythonEnv = new PythonEnv(PYTHON_VERSION, PYCUPPA_VENV_NAME);
    }

    private void removePycuppaTmpDir()
    {
        try {
            FileUtils.deleteDirectory(PYCUPPA_TMP_DIR);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to remove tmp dir: {}", PYCUPPA_TMP_DIR, e);
            System.exit(1);
        }
    }

    private void extractPycuppaToTmpDir()
    {
        try
        {
            if(PYCUPPA_TMP_DIR.exists())
            {
                CUP_LOGGER.debug("Removing existing tmp dir: {}", PYCUPPA_TMP_DIR);
                removePycuppaTmpDir();
            }

            if(PYCUPPA_RESOURCE_DIR.exists())
            {
                // This clause allows PycuppaExecutor to work when run from intellij
                CUP_LOGGER.debug("Copying {} from resource path({}) to tmp dir({})", PYCUPPA_PKG_NAME, PYCUPPA_RESOURCE_DIR, PYCUPPA_TMP_DIR);
                FileUtils.copyDirectory(PYCUPPA_RESOURCE_DIR, PYCUPPA_TMP_DIR);
            }
            else
            {
                String cuppaJarPath = new File( PythonEnv.class.getProtectionDomain().getCodeSource().getLocation().toURI() ).getPath();

                // Unzip with `jar` rather than `unzip` to avoid unzipping the whole jar
                String command = String.format("cd %s && jar -xf %s %s", TMP_DIR, cuppaJarPath, PYCUPPA_PKG_NAME);

                CUP_LOGGER.debug("Extracting {} from cuppa jar using command: {}", PYCUPPA_PKG_NAME, command);
                new BashCommand(command).logLevel(Level.DEBUG).showCommand().run();
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to extract {} to tmp dir: {}", PYCUPPA_PKG_NAME, PYCUPPA_TMP_DIR);
            System.exit(1);
        }
    }

    public void install()
    {
        try
        {
            mPythonEnv.initialize();

            if(!mPythonEnv.packageInstalled(PYCUPPA_PKG_NAME))
            {
                extractPycuppaToTmpDir();
                mPythonEnv.pipInstall(PYCUPPA_TMP_DIR.getPath(), true);
            }
            else
            {
                CUP_LOGGER.warn("Skipping installing {} as it already exists in virtual env({})", PYCUPPA_PKG_NAME, mPythonEnv.mVirtualEnvName);
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to install {} to virtual env({}): ", PYCUPPA_PKG_NAME, mPythonEnv.mVirtualEnvName);
            System.exit(1);
        }
        finally
        {
            CUP_LOGGER.debug("Removing tmp dir: " + PYCUPPA_TMP_DIR);
            removePycuppaTmpDir();
        }
    }

    public static void main(String[] args)
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
        PycuppaInstaller installer = new PycuppaInstaller();
        installer.install();
    }
}
