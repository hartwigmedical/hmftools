package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.Nullable;

public abstract class ShellCommand
{
    public final ProcessBuilder mProcessBuilder;

    @Nullable
    private Level mLogLevel = Level.INFO;
    private boolean mShowCommand = false;

    private int mExitCode = 0;
    private List<String> mStdout;

    public ShellCommand(ProcessBuilder processBuilder)
    {
        mProcessBuilder = processBuilder;
    }

    public ShellCommand logLevel(@Nullable Level logLevel)
    {
        mLogLevel = logLevel;
        return this;
    }

    public ShellCommand showCommand()
    {
        mShowCommand = true;
        return this;
    }

    public abstract String toString();

    public ShellCommand run()
    {
        try
        {
            if(mShowCommand)
            {
                CUP_LOGGER.info("Running command: " + this);
            }

            Process process = mProcessBuilder
                    .redirectErrorStream(true) // Merge stdout and stderr
                    .start();

            mStdout = new ArrayList<>();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String line;
            while( (line = reader.readLine()) != null )
            {
                mStdout.add(line);
                if(mLogLevel != null)
                {
                    CUP_LOGGER.log(mLogLevel, line);
                }
            }

            mExitCode = process.waitFor();
            if(mExitCode > 0)
            {
                throw new RuntimeException();
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("command({}) failed with exit code {}: {}", mProcessBuilder.command(), String.valueOf(mExitCode), e);
            System.exit(mExitCode);
        }

        return this;
    }

    public List<String> getStdout(){ return mStdout; }
}
