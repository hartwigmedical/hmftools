package com.hartwig.hmftools.common.progress;

import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

import javax.annotation.concurrent.ThreadSafe;

import com.google.common.base.Strings;
import com.google.common.util.concurrent.AtomicDouble;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@ThreadSafe
public class FutureProgressTracker {

    private static final Logger LOGGER = LogManager.getLogger(FutureProgressTracker.class);

    private final double increment;
    private final int minExpected;
    private final AtomicInteger expected = new AtomicInteger();
    private final AtomicInteger complete = new AtomicInteger();
    private final AtomicDouble previousPercentComplete = new AtomicDouble();

    public FutureProgressTracker() {
       this(0.1, 0);
    }

    public FutureProgressTracker(final double increment, final int minExpected) {
        this.increment = increment;
        this.minExpected = minExpected;
    }

    public Runnable add(@NotNull final Runnable callable) {
        expected.incrementAndGet();

        return () -> {
            callable.run();
            completed();
        };
    }

    public <T> Callable<T> add(@NotNull final Callable<T> callable) {
        expected.incrementAndGet();

        return () -> {
            T result = callable.call();
            completed();
            return result;
        };
    }

    private void completed() {
        double completeInt = complete.incrementAndGet();
        double expectedInt = expected.get();
        if (expectedInt >= minExpected) {
            double percentComplete = completeInt / expectedInt;
            synchronized (previousPercentComplete) {
                if (completeInt == expectedInt || percentComplete > previousPercentComplete.get() + increment) {
                    LOGGER.info("{}", complete(percentComplete));
                    previousPercentComplete.set(percentComplete);
                }
            }
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

