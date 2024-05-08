package com.hartwig.hmftools.cup.cli;

import java.util.StringJoiner;

public class BashCommand extends ShellCommand
{
    private final String mCommand;

    public BashCommand(String command)
    {
        super(new ProcessBuilder("bash", "-c", command));
        mCommand = command;
    }

    private static String formCommand(String[] args)
    {
        StringJoiner joiner = new StringJoiner(" ");

        for(String arg : args)
            joiner.add(arg);

        return joiner.toString();
    }

    public BashCommand(String... args)
    {
        this(formCommand(args));
    }

    @Override
    public String toString()
    {
        return String.format("bash -c \"%s\"", mCommand);
    }
}

