package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import com.google.common.base.Strings;

import org.jetbrains.annotations.NotNull;

public class AmberTaskCompletion
{
    private final int mExpected;
    private double mPreviousPercentComplete;

    public AmberTaskCompletion(int numTasks)
    {
        mExpected = numTasks;
        mPreviousPercentComplete = 0;
    }

    public void progress(int numRemainingTasks)
    {
        double percentComplete = 1.0 - (double)numRemainingTasks / mExpected;
        if (percentComplete - mPreviousPercentComplete >= 0.05 ||
            mPreviousPercentComplete != 1.0 && percentComplete == 1.0)
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