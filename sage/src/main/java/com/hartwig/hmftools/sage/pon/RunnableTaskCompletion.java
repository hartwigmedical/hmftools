package com.hartwig.hmftools.sage.pon;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import com.google.common.base.Strings;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class RunnableTaskCompletion
{
    private int expected = 0;
    private int complete = 0;
    private double previousPercentComplete = 0;

    @NotNull
    Runnable task(@NotNull final Runnable runnable)
    {
        expected++;
        return () ->
        {
            runnable.run();
            completed();
        };
    }

    private void completed()
    {
        double percentComplete = ((double) ++complete) / expected;
        if(expected == complete || percentComplete > previousPercentComplete + 0.1)
        {
            SG_LOGGER.info("{}", complete(percentComplete));
            previousPercentComplete = percentComplete;
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
