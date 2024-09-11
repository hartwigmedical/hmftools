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

    @Nullable
    private long mTimeout = Long.MAX_VALUE;
    private TimeUnit mTimeoutUnit = TimeUnit.SECONDS;

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

            mComplete = process.waitFor(mTimeout, mTimeoutUnit);
            CUP_LOGGER.info("Command complete: [{}]", this);
            mExitCode = process.exitValue();
        }
        catch(InterruptedException | IllegalThreadStateException e)
        {
            CUP_LOGGER.error("Failed to run command({}) due to time out after {} {}", mProcessBuilder.command(), mTimeout, mTimeoutUnit.toString().toLowerCase());
            mExitCode = 1;
        }
        catch(IOException | RuntimeException e)
        {
            CUP_LOGGER.error("Failed to run command({})", mProcessBuilder.command());
            mExitCode = 1;
        }

        return this;
    }

    @Override
    public String toString()
    {
        return String.join(" ", mProcessBuilder.command());
    }
}
