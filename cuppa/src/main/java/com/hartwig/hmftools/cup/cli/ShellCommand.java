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
    private Level mLogLevel = Level.INFO;
    private boolean mShowCommand = true;
    private long mTimeout = Long.MAX_VALUE;
    private TimeUnit mTimeoutUnit = TimeUnit.SECONDS;

    private int mExitCode;
    private List<String> mStdout;
    private boolean mComplete = false;

    private static final String THREAD_NAME = "bash";

    public ShellCommand(ProcessBuilder processBuilder)
    {
        mProcessBuilder = processBuilder;
    }

    private void checkComplete()
    {
        if(mComplete)
            throw new IllegalStateException("Cannot set options or run a command that is already complete");
    }

    public ShellCommand logLevel(@Nullable Level logLevel)
    {
        checkComplete();
        mLogLevel = logLevel;
        return this;
    }

    public ShellCommand showCommand()
    {
        checkComplete();
        mShowCommand = true;
        return this;
    }

    public ShellCommand timeout(long timeout, TimeUnit unit)
    {
        checkComplete();
        mTimeout = timeout;
        mTimeoutUnit = unit;
        return this;
    }

    public ShellCommand timeout(long timeout)
    {
        checkComplete();
        mTimeout = timeout;
        return this;
    }

    public ShellCommand run()
    {
        checkComplete();
        try
        {
            if(mShowCommand)
                CUP_LOGGER.info("Running command: [{}]", this);

            Process process = mProcessBuilder
                    .redirectErrorStream(true) // Merge stdout and stderr
                    .start();

            mComplete = process.waitFor(mTimeout, mTimeoutUnit);
            CUP_LOGGER.info("Command [{}] complete", this);
            mStdout = captureStdout(process, mLogLevel);
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

    private static List<String> captureStdout(Process process, Level logLevel) {
        List<String> stdout = new ArrayList<>();
        try {
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));

            String line;
            while ((line = reader.readLine()) != null) {
                stdout.add(line);
                if (logLevel != null)
                    CUP_LOGGER.log(logLevel, line);
            }
            return stdout;
        } catch (IOException e) {
            throw new RuntimeException("Failed to caputre output from process!", e);
        }
    }

    public List<String> getStdout(boolean warnEmpty)
    {
        if(warnEmpty & mStdout.isEmpty())
            CUP_LOGGER.warn("stdout is empty for command({})", mProcessBuilder.command());

        return mStdout;
    }

    public List<String> getStdout()
    {
        return getStdout(true);
    }

    @Override
    public String toString()
    {
        return String.join(" ", mProcessBuilder.command());
    }
}
