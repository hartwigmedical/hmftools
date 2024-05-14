package com.hartwig.hmftools.cup.cli;

public class PythonEnvCommand extends ShellCommand
{
    private final String mCommand;

    public PythonEnvCommand(PythonEnv pythonEnvironment, String command)
    {
        super(new ProcessBuilder("bash", "-c", formCommand(pythonEnvironment, command)));

        mCommand = formCommand(pythonEnvironment, command);
    }

    private static String formCommand(PythonEnv pythonEnvironment, String command)
    {
        return String.format("%s && source %s/bin/activate && %s && deactivate",
                pythonEnvironment.exportPyenvRootCommand(),
                pythonEnvironment.virtualEnvPath(),
                command
        );
    }

    @Override
    public String toString()
    {
        return mCommand;
    }
}
