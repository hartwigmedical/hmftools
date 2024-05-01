package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

public class PythonEnvCommand extends ShellCommand
{
    public PythonEnvCommand(PythonEnv pythonEnvironment, String command)
    {
        super(createProcessBuilder(pythonEnvironment, command));
    }

    public static ProcessBuilder createProcessBuilder(PythonEnv pythonEnvironment, String command)
    {
        String bashCommand =  String.format("source %s/bin/activate && %s && deactivate", pythonEnvironment.virtualEnvPath(), command);
        return new ProcessBuilder("bash", "-c", bashCommand);
    }

    @Override
    public String toString()
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
        return mProcessBuilder.command().toString();
    }
}
