package com.hartwig.hmftools.cup.cli;

import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class PythonInterpreter
{
    public final File mBinaryPath;

    public PythonInterpreter(final String pythonBinaryPath)
    {
        checkValid(pythonBinaryPath);
        mBinaryPath = new File(pythonBinaryPath);
    }

    private void checkValid(final String pythonBinaryPath)
    {
        Path path = Path.of(pythonBinaryPath);
        if(!Files.isExecutable(path) || !path.getFileName().startsWith("python"))
        {
            CUP_LOGGER.error("Invalid python interpreter: {}", pythonBinaryPath);
            System.exit(1);
        }
    }

    private class PythonCommand extends ShellCommand
    {
        public PythonCommand(final String command)
        {
            super(new ProcessBuilder("bash", "-c", mBinaryPath + " " + command));
        }
    }

    public PythonCommand command(final String command)
    {
        return new PythonCommand(command);
    }

    public PythonCommand command(final String... args)
    {
        return new PythonCommand(String.join(" ", args));
    }

    public String getVersion()
    {
        return command("-c 'import platform; print(platform.python_version())'")
            .logLevel(null).timeout(10).run()
            .getStdout().get(0);
    }

    private List<String> getInstalledPackages(boolean editableOnly)
    {
        List<String> installedPackages = new ArrayList<>();

        try
        {
            String commandString = "-m pip list --format=freeze --disable-pip-version-check";
            commandString += (editableOnly) ? " --editable" : " --exclude-editable";

            List<String> stdout = command(commandString)
                    .logLevel(null).timeout(10)
                    .run().getStdout(false);

            for(String line : stdout)
            {
                String packageName = line.split("==")[0];
                installedPackages.add(packageName);
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to get installed python packages: ", e);
        }

        return installedPackages;
    }

    public PythonInterpreter requirePackages(String... packageNames)
    {
        List<String> nonEditablePackages = getInstalledPackages(false);
        List<String> editablePackages = getInstalledPackages(true);

        int missingNonEditiblePackages = 0;
        for(String packageName : packageNames)
        {
            if(editablePackages.contains(packageName))
            {
                CUP_LOGGER.error("Python package({}) is installed in editable mode and does not work when called with PythonInterpreter", packageName);
                missingNonEditiblePackages++;
                continue;
            }

            if(!nonEditablePackages.contains(packageName))
            {
                CUP_LOGGER.error("Python package({}) missing", packageName);
                missingNonEditiblePackages++;
            }
        }

        if(missingNonEditiblePackages > 0)
            System.exit(1);

        return this;
    }
}
