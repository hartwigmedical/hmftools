package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.File;
import java.util.List;

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

    public final PythonEnvironment mPythonEnvironment;

    public PycuppaInstaller(String virtualEnvPath)
    {
        mPythonEnvironment = new PythonEnvironment(PYTHON_VERSION, virtualEnvPath);
    }

    private void removePycuppaTmpDir()
    {
        try {
            FileUtils.deleteDirectory(PYCUPPA_TMP_DIR);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to remove tmp dir({}): ", PYCUPPA_TMP_DIR, e);
            System.exit(1);
        }
    }

    private void extractPycuppaToTmpDir()
    {
        try
        {
            if(PYCUPPA_TMP_DIR.exists())
            {
                CUP_LOGGER.debug("Removing existing pycuppa tmp dir: " + PYCUPPA_TMP_DIR);
                removePycuppaTmpDir();
            }

            if(PYCUPPA_RESOURCE_DIR.exists())
            {
                // This clause allows PycuppaExecutor to work when run from intellij
                CUP_LOGGER.debug("Copying pycuppa from resource path({}) to tmp dir({})", PYCUPPA_RESOURCE_DIR, PYCUPPA_TMP_DIR);
                FileUtils.copyDirectory(PYCUPPA_RESOURCE_DIR, PYCUPPA_TMP_DIR);
            }
            else
            {
                String cuppaJarPath = new File( PythonEnvironment.class.getProtectionDomain().getCodeSource().getLocation().toURI() ).getPath();

                // Unzip with `jar` rather than `unzip` to avoid unzipping the whole jar
                String command = String.format("cd %s && jar -xf %s pycuppa", TMP_DIR, cuppaJarPath);

                CUP_LOGGER.debug("Extracting pycuppa from cuppa jar using command: " + command);
                new BashCommand(command).logLevel(Level.DEBUG).showCommand().run();
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to extract pycuppa to tmp dir({})", PYCUPPA_TMP_DIR, e);
            System.exit(1);
        }
    }

    private boolean pycuppaInstalled()
    {
        boolean pycuppaExists = false;

        List<String> stdout = new PythonCommand(mPythonEnvironment, "python -m pip --disable-pip-version-check list")
                .logLevel(null).run().getStdout();

        for(String line : stdout)
        {
            if(line.startsWith("pycuppa"))
            {
                pycuppaExists = true;
                break;
            }
        }

        return pycuppaExists;
    }

    public void install()
    {
        try
        {
            mPythonEnvironment.initialize(true);

            if(!pycuppaInstalled())
            {
                extractPycuppaToTmpDir();
                mPythonEnvironment.pipUpgrade();
                mPythonEnvironment.pipInstall(PYCUPPA_TMP_DIR.getPath());
            }
            else
            {
                CUP_LOGGER.warn("Skipping installing pycuppa as it already exists in the virtual environment: " + mPythonEnvironment.mVirtualEnvPath);
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to install pycuppa to virtual env({}): ", mPythonEnvironment.mVirtualEnvPath, e);
            System.exit(1);
        }
        finally
        {
            CUP_LOGGER.debug("Removing pycuppa tmp dir: " + PYCUPPA_TMP_DIR);
            removePycuppaTmpDir();
        }
    }

    public static void main(String[] args)
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
        PycuppaInstaller installer = new PycuppaInstaller(args[0]);
        installer.install();
    }
}
