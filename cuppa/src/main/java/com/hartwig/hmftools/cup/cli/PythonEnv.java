package com.hartwig.hmftools.cup.cli;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;

public class PythonEnv
{
    public static final String DEFAULT_PYTHON_VERSION = "3.9.4";
    public static final String DEFAULT_PYENV_DIR = System.getProperty("user.home") + "/.pyenv";

    public final String mPythonVersion;
    public final File mVirtualEnvName;
    public final File mPyenvDir;

    public PythonEnv(String pythonVersion, String virtualEnvName, String pyenvDir)
    {
        mVirtualEnvName = new File(virtualEnvName);
        mPythonVersion = pythonVersion;
        mPyenvDir = new File(pyenvDir);
    }

    public PythonEnv(String pythonVersion, String virtualEnvName)
    {
        this(pythonVersion, virtualEnvName, DEFAULT_PYENV_DIR);
    }

    public File pyenvPath(){ return new File(mPyenvDir + "/libexec/pyenv"); }

    public File virtualEnvPath() { return new File(String.format("%s/versions/%s/envs/%s", mPyenvDir, mPythonVersion, mVirtualEnvName)); }

    public File pythonPath() { return new File(virtualEnvPath() + "/bin/python"); }

    public String exportPyenvRootCommand(){ return "export PYENV_ROOT="+mPyenvDir; }

    @VisibleForTesting
    void installPyenv()
    {
        if(!mPyenvDir.exists())
        {
            // This command also installs pyenv-virtualenv to $PYENV_ROOT/plugins/pyenv-virtualenv
            ShellCommand command = new BashCommand(exportPyenvRootCommand(), "&& curl https://pyenv.run | bash")
                    .logLevel(Level.DEBUG);

            CUP_LOGGER.info("Installing pyenv with command: " + command);
            command.run();
        }
        else
        {
            CUP_LOGGER.warn("Skipping installing pyenv as it exists at: " + mPyenvDir);
        }

        if(!mPyenvDir.exists())
        {
            CUP_LOGGER.error("Failed to install pyenv to " + mPyenvDir);
            System.exit(1);
        }
    }

    @VisibleForTesting
    File getRcFilePath() throws FileNotFoundException
    {
        String[] basenamesToTry = {".zshrc", ".bashrc"};

        String pathToRcFile = null;
        for(String basename : basenamesToTry)
        {
            String path = System.getProperty("user.home") + "/" + basename;
            if(new File(path).isFile())
                pathToRcFile = path;
        }

        if(pathToRcFile == null)
            throw new FileNotFoundException("Failed to locate one of the following rc files: " + String.join(", ", basenamesToTry));

        return new File(pathToRcFile);
    }

    public void addPyenvPathsToRcFile()
    {
        try
        {
            File RcFilePath = getRcFilePath();

            String lines = String.join("\n",
                "# Pyenv paths",
                "export PYENV_ROOT="+mPyenvDir,
                "[[ -d $PYENV_ROOT/bin ]] && export PATH=\"$PYENV_ROOT/bin:$PATH\"",
                "eval \"$(pyenv init -)\"",
                "eval \"$(pyenv virtualenv-init -)\""
            );

            boolean linesExist = FileUtils.readFileToString(RcFilePath, "UTF-8").contains(lines);
            if(!linesExist)
            {
                CUP_LOGGER.info("Adding pyenv paths to rc file: " + RcFilePath);
                String command = String.format("echo -e '\n%s' >> %s", lines, RcFilePath);
                new BashCommand(command).logLevel(null).run();
            }
            else
            {
                CUP_LOGGER.warn("Skipping adding pyenv paths to rc file: " + RcFilePath);
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.warn("Failed to add pyenv paths to rc file (using pyenv in interactive shell will not work): " + e);
        }
    }

    @VisibleForTesting
    void installPython()
    {
        File pythonPathInPyenv = new File(String.format("%s/versions/%s/bin/python", mPyenvDir, mPythonVersion));

        if(!pythonPathInPyenv.exists())
        {
            ShellCommand command = new BashCommand(exportPyenvRootCommand(), "&&", pyenvPath().toString(), "install", mPythonVersion)
                    .logLevel(Level.DEBUG);

            CUP_LOGGER.info("Installing python version {} with command: {}", mPythonVersion, command);
            command.run();
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
            ShellCommand command = new BashCommand(
                    exportPyenvRootCommand(), "&&",
                    pyenvPath().toString(), "virtualenv", mPythonVersion, mVirtualEnvName.toString()
            ).logLevel(Level.DEBUG);

            CUP_LOGGER.info("Creating python virtual environment with command: " + command);
            command.run();
        }
        else
        {
            CUP_LOGGER.warn("Skipping creating python virtual environment as binary exists at: " + pythonPath());
        }
    }

    @VisibleForTesting
    void removeExistingPyenv()
    {
        CUP_LOGGER.info("Removing existing pyenv at: " + mPyenvDir);
        try
        {
            FileUtils.deleteDirectory(mPyenvDir);
        }
        catch(IOException e)
        {
            CUP_LOGGER.warn("Failed to remove pyenv: " + e);
        }
    }

    public PythonEnv initialize()
    {
        if(pythonPath().exists())
            return this;

        installPyenv();
        installPython();
        createVirtualEnvironment();
        return this;
    }

    private String pipList()
    {
        return new PythonEnvCommand(this, "pip --disable-pip-version-check list")
                .logLevel(null).run().getStdoutAsString();
    }

    public boolean packageInstalled(String packageName)
    {
        String stdout = pipList();
        return stdout.contains(packageName);
    }

    public void pipInstall(String args, boolean upgradePip)
    {
        String command = "pip install " + args;
        if(upgradePip)
            command = "pip install --upgrade pip && " + command;

        CUP_LOGGER.info("Installing packages with command: " + command);
        new PythonEnvCommand(this, command).logLevel(Level.DEBUG).run();
    }

    public void pipInstall(String args){ pipInstall(args, true); }

    public PythonEnv checkRequiredPackages(String[] packageNames)
    {
        if(!pythonPath().exists())
        {
            CUP_LOGGER.error("Missing or invalid python virtual environment: " + virtualEnvPath());
            System.exit(1);
        }

        String stdout = pipList();

        List<String> missingPackages = new ArrayList<>();
        for(String packageName : packageNames)
        {
            if(!stdout.contains(packageName))
                missingPackages.add(packageName);
        }

        if(missingPackages.size() > 0)
        {
            CUP_LOGGER.error("Python virtual environment({}) missing the following packages: {}",
                    virtualEnvPath(), String.join(", ", missingPackages));
            System.exit(1);
        }

        return this;
    }

    public PythonEnv checkRequiredPackage(String packageName)
    {
        return checkRequiredPackages(new String[]{packageName});
    }
}
