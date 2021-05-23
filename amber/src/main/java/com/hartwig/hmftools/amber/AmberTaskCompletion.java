package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.concurrent.Callable;

import com.google.common.base.Strings;

import org.jetbrains.annotations.NotNull;

public class AmberTaskCompletion
{
    private int mExpected;
    private int mComplete;
    private double mPreviousPercentComplete;

    public AmberTaskCompletion()
    {
        mExpected = 0;
        mComplete = 0;
        mPreviousPercentComplete = 0;
    }

    public <T> Callable<T> task(@NotNull final Callable<T> callable)
    {
        mExpected++;

        return () ->
        {
            T result = callable.call();
            completed();
            return result;
        };
    }

    private void completed()
    {
        double percentComplete = ((double) ++mComplete) / mExpected;
        if(mExpected == mComplete || percentComplete > mPreviousPercentComplete + 0.1)
        {
            AMB_LOGGER.info("{}", complete(percentComplete));
            mPreviousPercentComplete = percentComplete;
        }
    }

    @NotNull
    private static String complete(double percent)
    {
        int roundedPercent = (int) Math.round(percent * 100);
        int hashCount = Math.min(20, roundedPercent / 5);
        int gapCount = Math.max(0, 20 - hashCount);

        return "  [" + Strings.repeat("#", hashCount) + Strings.repeat(" ", gapCount) + "] " + roundedPercent + "% complete";
    }
}
