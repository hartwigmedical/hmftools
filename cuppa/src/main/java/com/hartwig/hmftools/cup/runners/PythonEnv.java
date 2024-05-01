package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.File;
import java.io.IOException;

import com.google.common.annotations.VisibleForTesting;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;

public class PythonEnv
{
    public static final String DEFAULT_PYTHON_VERSION = "3.9.4";

    private static final File PYENV_DIR = new File(System.getProperty("user.home") + "/.pyenv");
    private static final File PYENV_PATH = new File(PYENV_DIR + "/libexec/pyenv");

    public final String mPythonVersion;
    public final File mVirtualEnvName;

    public PythonEnv(String pythonVersion, String virtualEnvName)
    {
        mVirtualEnvName = new File(virtualEnvName);
        mPythonVersion = pythonVersion;
    }

    public File virtualEnvPath() { return new File(String.format("%s/versions/%s/envs/%s", PYENV_DIR, mPythonVersion, mVirtualEnvName)); }

    public File pythonPath() { return new File(virtualEnvPath() + "/bin/python"); }

    @VisibleForTesting
    void installPyenv()
    {
        if(!PYENV_DIR.exists())
        {
            String command;

            command = "curl https://pyenv.run | bash"; // This also installs pyenv-virtualenv to $HOME/.pyenv/plugins/pyenv-virtualenv
            CUP_LOGGER.info("Installing pyenv with command: " + command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else
        {
            CUP_LOGGER.warn("Skipping installing pyenv as it exists at: " + PYENV_DIR);
        }

        if(!PYENV_DIR.exists())
        {
            CUP_LOGGER.error("Failed to install pyenv to " + PYENV_DIR);
            System.exit(1);
        }
    }

    @VisibleForTesting
    void addPyenvPathsToRcFile()
    {
        String shellPath = new BashCommand("echo $SHELL").logLevel(null).run().getStdout().get(0);
        String shellName = new File(shellPath).getName();
        File RcFilePath = new File(String.format("%s/.%src", System.getProperty("user.home"), shellName));

        if(!RcFilePath.isFile())
        {
            CUP_LOGGER.error("Fail to find rc file({}) to add pyenv paths. Using pyenv in interactive shell will not work", RcFilePath);
            System.exit(1);
        }

        try
        {
            String lines = String.join("\n",
                "# Pyenv paths",
                "export PYENV_ROOT="+ PYENV_DIR,
                "[[ -d $PYENV_ROOT/bin ]] && export PATH=\"$PYENV_ROOT/bin:$PATH\"",
                "eval \"$(pyenv init -)\"",
                "eval \"$(pyenv virtualenv-init -)\""
            );

            boolean linesExist = FileUtils.readFileToString(RcFilePath, "UTF-8").contains(lines);
            if(!linesExist)
            {
                String command;

                CUP_LOGGER.info("Adding pyenv paths to rc file: " + RcFilePath);
                command = String.format("echo -e '\n%s' >> %s", lines, RcFilePath);
                new BashCommand(command).logLevel(null).run();
            }
            else
            {
                CUP_LOGGER.warn("Skipping adding pyenv paths to rc file: " + RcFilePath);
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to update rc file({}): {}", RcFilePath, e);
            System.exit(1);
        }
    }

    @VisibleForTesting
    void installPython()
    {
        File pythonPathInPyenv = new File(String.format("%s/versions/%s/bin/python", PYENV_DIR, mPythonVersion));

        if(!pythonPathInPyenv.exists())
        {
            String command = String.format("%s install %s", PYENV_PATH, mPythonVersion);
            CUP_LOGGER.info("Installing python version {} with command: {}", mPythonVersion, command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else
        {
            CUP_LOGGER.warn("Skipping installing python {} as it exists at: {}", mPythonVersion, pythonPathInPyenv);
        }
    }

    @VisibleForTesting
    void createVirtualEnvironment()
    {
        if(!pythonPath().exists())
        {
            String command = String.format("%s virtualenv %s %s", PYENV_PATH, mPythonVersion, mVirtualEnvName);
            CUP_LOGGER.info("Creating python virtual environment with command: " + command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else
        {
            CUP_LOGGER.warn("Skipping creating python virtual environment as binary exists at: " + pythonPath());
        }
    }

    @VisibleForTesting
    void removeExistingPyenv()
    {
        CUP_LOGGER.info("Removing existing pyenv at: " + PYENV_DIR);
        try {
            FileUtils.deleteDirectory(PYENV_DIR);
        } catch(IOException e) {
            CUP_LOGGER.warn("Failed to remove pyenv: " + e);
        }
    }

    public PythonEnv initialize()
    {
        if(pythonPath().exists())
            return this;

        installPyenv();
        addPyenvPathsToRcFile();
        installPython();
        createVirtualEnvironment();
        return this;
    }

    public void pipUpgrade()
    {
        String command = "python -m pip install --upgrade pip";
        CUP_LOGGER.info("Upgrading pip with command: " + command);
        new PythonEnvCommand(this, command).showCommand().logLevel(Level.DEBUG).run();
    }

    public void pipInstall(String args)
    {
        String command = "python -m pip install " + args;
        CUP_LOGGER.info("Installing packages with command: " + command);
        new PythonEnvCommand(this, command).showCommand().logLevel(Level.DEBUG).run();
    }
}
