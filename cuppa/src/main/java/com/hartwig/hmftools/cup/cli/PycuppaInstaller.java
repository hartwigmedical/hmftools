package com.hartwig.hmftools.cup.cli;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.LOG_LEVEL;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.LOG_LEVEL_DESC;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.cli.PythonEnv.DEFAULT_PYENV_DIR;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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

    public static final String INSTALL_DIR = "install_dir";
    public static final String INSTALL_DIR_DESC = "Path to the pyenv dir in which to install pycuppa (will be created if not existing). Default: " + DEFAULT_PYENV_DIR;

    public static final String UPDATE_RC_FILE = "update_rc_file";
    public static final String UPDATE_RC_FILE_DESC = "Add pyenv paths and initialization logic to bashrc or zshrc file? Only required for use of pycuppa in an interactive shell";

    public final ConfigBuilder mConfig;
    public final PythonEnv mPythonEnv;

    public PycuppaInstaller(ConfigBuilder config)
    {
        mConfig = config;

        Configurator.setLevel(CUP_LOGGER.getName(), Level.valueOf(config.getValue(LOG_LEVEL)));

        String installDir = getInstallDir();
        mPythonEnv = new PythonEnv(PYTHON_VERSION, PYCUPPA_VENV_NAME, installDir, false);
    }

    private String getInstallDir()
    {
        String installDir = null;
        if(!mConfig.hasValue(INSTALL_DIR))
        {
            CUP_LOGGER.info("Installing to default pyenv dir: {}", DEFAULT_PYENV_DIR);
            installDir = DEFAULT_PYENV_DIR;
        }
        else
        {
            String rawInstallDir = mConfig.getValue(INSTALL_DIR);
            try
            {
                installDir = new File(rawInstallDir).getCanonicalPath();
            }
            catch(IOException e)
            {
                CUP_LOGGER.error("Failed to get absolute path to installDir({}): {}", rawInstallDir, e.toString());
                System.exit(1);
            }
        }

        return installDir;
    }

    private void removePycuppaTmpDir()
    {
        try
        {
            FileUtils.deleteDirectory(PYCUPPA_TMP_DIR);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to remove tmp dir: {}", PYCUPPA_TMP_DIR, e);
            System.exit(1);
        }
    }

    private void extractPycuppa()
    {
        try
        {
            if(PYCUPPA_TMP_DIR.exists())
            {
                CUP_LOGGER.debug("Removing existing tmp dir: {}", PYCUPPA_TMP_DIR);
                removePycuppaTmpDir();
            }

            String cuppaJarPath = new File( PythonEnv.class.getProtectionDomain().getCodeSource().getLocation().toURI() ).getPath();

            // Unzip with `jar` rather than `unzip` to avoid unzipping the whole jar
            ShellCommand command = new BashCommand(String.format("cd %s && jar -xf %s %s", TMP_DIR, cuppaJarPath, PYCUPPA_PKG_NAME))
                    .logLevel(Level.DEBUG).showCommand();

            CUP_LOGGER.debug("Extracting {} from cuppa jar using command: {}", PYCUPPA_PKG_NAME, command);
            command.run();
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to extract {} to tmp dir({}): {}", PYCUPPA_PKG_NAME, PYCUPPA_TMP_DIR, e.toString());
            System.exit(1);
        }
    }

    public void install()
    {
        try
        {
            mPythonEnv.install();

            if(mConfig.hasFlag(UPDATE_RC_FILE))
                mPythonEnv.updateRcFile();

            if(mPythonEnv.packageInstalled(PYCUPPA_PKG_NAME))
            {
                CUP_LOGGER.warn("Skipping installing {} as it already exists in virtual env({})", PYCUPPA_PKG_NAME, mPythonEnv.mVirtualEnvName);
                return;
            }

            File pycuppaDir;
            if(PYCUPPA_RESOURCE_DIR.exists())
            {
                // This clause allows pycuppa to be installed when run from intellij
                CUP_LOGGER.debug("Using pycuppa from resource path: {}", PYCUPPA_RESOURCE_DIR);
                pycuppaDir = PYCUPPA_RESOURCE_DIR;
            }
            else
            {
                extractPycuppa();
                pycuppaDir = PYCUPPA_TMP_DIR;
            }
            mPythonEnv.pipInstall(true, pycuppaDir.toString());

            CUP_LOGGER.info("Completed installation of " + PYCUPPA_PKG_NAME);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to install {} to virtual env({}): {}", PYCUPPA_PKG_NAME, mPythonEnv.mVirtualEnvName, e);
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
        ConfigBuilder config = new ConfigBuilder(APP_NAME);
        config.addConfigItem(INSTALL_DIR, false, INSTALL_DIR_DESC);
        config.addFlag(UPDATE_RC_FILE, UPDATE_RC_FILE_DESC);
        config.addConfigItem(LOG_LEVEL, false, LOG_LEVEL_DESC, Level.DEBUG.toString());
        config.checkAndParseCommandLine(args);

        PycuppaInstaller installer = new PycuppaInstaller(config);
        installer.install();
    }
}