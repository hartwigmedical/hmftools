package com.hartwig.hmftools.cdr3;

import com.google.common.base.Strings;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BamTaskCompletion
{
    private static final Logger logger = LogManager.getLogger(BamTaskCompletion.class);

    private final int mExpected;
    private double mPreviousPercentComplete;

    public BamTaskCompletion(int numTasks)
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
            logger.info("{}", complete(percentComplete));
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