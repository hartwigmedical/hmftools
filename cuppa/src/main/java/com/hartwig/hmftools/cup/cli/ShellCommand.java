package com.hartwig.hmftools.cup.cli;

import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.Nullable;

public abstract class ShellCommand
{
    public final ProcessBuilder mProcessBuilder;

    private int mExitCode;
    private boolean mComplete = false;

    public ShellCommand(ProcessBuilder processBuilder)
    {
        mProcessBuilder = processBuilder;
    }

    private void checkComplete()
    {
        if(mComplete)
            throw new IllegalStateException("Cannot set options or run a command that is already complete");
    }

    public ShellCommand run()
    {
        checkComplete();

        try
        {
            CUP_LOGGER.info("Running command: [{}]", this);

            Process process = mProcessBuilder
                    .redirectErrorStream(true) // Merge stdout and stderr
                    .inheritIO()
                    .start();

            mComplete = process.waitFor(Long.MAX_VALUE, TimeUnit.SECONDS);
            mExitCode = process.exitValue();

            if(mExitCode != 0)
            {
                CUP_LOGGER.error("Failed to run command: {}", mProcessBuilder.command());
                throw new RuntimeException();
            }

            CUP_LOGGER.info("Command complete: [{}]", this);
        }
        catch(IOException | InterruptedException e)
        {
            CUP_LOGGER.error("Failed to run command: {}", mProcessBuilder.command());
            throw new RuntimeException();
        }

        return this;
    }

    @Override
    public String toString()
    {
        return String.join(" ", mProcessBuilder.command());
    }
}
