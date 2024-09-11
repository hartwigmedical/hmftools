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
            super(new ProcessBuilder("/bin/bash", "-c", mBinaryPath + " " + command));
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
}
