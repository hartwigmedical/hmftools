package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

public class PythonCommand extends ShellCommand
{
    public PythonCommand(PythonEnvironment pythonEnvironment, String command)
    {
        super(createProcessBuilder(pythonEnvironment, command));
    }

    public static ProcessBuilder createProcessBuilder(PythonEnvironment pythonEnvironment, String command)
    {
        if(!command.startsWith("python"))
        {
            throw new IllegalStateException("Invalid python command: " + command);
        }

        String bashCommand =  String.format("source %s/bin/activate && %s && deactivate", pythonEnvironment.mVirtualEnvPath, command);
        return new ProcessBuilder("bash", "-c", bashCommand);
    }

    @Override
    public String toString()
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
        return mProcessBuilder.command().toString();
    }
}
