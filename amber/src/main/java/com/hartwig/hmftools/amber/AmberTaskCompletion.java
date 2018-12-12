package com.hartwig.hmftools.amber;

import java.util.concurrent.Callable;

import com.google.common.base.Strings;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class AmberTaskCompletion {

    private static final Logger LOGGER = LogManager.getLogger(AmberApplication.class);

    private int expected = 0;
    private int complete = 0;
    private double previousPercentComplete = 0;

    public <T> Callable<T> task(@NotNull final Callable<T> callable) {
        expected++;

        return () -> {
            T result = callable.call();
            completed();
            return result;
        };
    }

    private void completed() {
        double percentComplete = ((double) ++complete) / expected;
        if (expected == complete || percentComplete > previousPercentComplete + 0.1) {
            LOGGER.info("{}", complete(percentComplete));
            previousPercentComplete = percentComplete;
        }

    }

    @NotNull
    private static String complete(double percent) {
        int roundedPercent = (int) Math.round(percent * 100);
        int hashCount = Math.min(20, roundedPercent / 5);
        int gapCount = Math.max(0, 20 - hashCount);

        return "  [" + Strings.repeat("#", hashCount) + Strings.repeat(" ", gapCount) + "] " + roundedPercent + "% complete";
    }
}
