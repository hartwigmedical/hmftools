package com.hartwig.hmftools.cup.runners;

public class BashCommand extends ShellCommand
{
    private final String mCommand;

    public BashCommand(String command)
    {
        super(new ProcessBuilder("bash", "-c", command));
        mCommand = command;
    }

    @Override
    public String toString()
    {
        return String.format("bash -c \"%s\"", mCommand);
    }
}

