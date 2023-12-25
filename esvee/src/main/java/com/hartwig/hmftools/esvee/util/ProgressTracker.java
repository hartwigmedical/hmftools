package com.hartwig.hmftools.esvee.util;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class ProgressTracker implements AutoCloseable
{
    private final ScheduledExecutorService mReporter;

    @Nullable
    private final String mDescription;
    private final long mStartTimeNanos;

    private final int mTaskCount;
    private final AtomicInteger mTasksComplete = new AtomicInteger();
    private final Deque<Long> mPreviousTimeNanos = new ArrayDeque<>();
    private final Deque<Integer> mPreviousTasks = new ArrayDeque<>();

    public ProgressTracker(final int maximum, final int reportFrequencySeconds)
    {
        this(null, maximum, reportFrequencySeconds);
    }

    public ProgressTracker(@Nullable final String description, final int maximum, final int reportFrequencySeconds)
    {
        mDescription = description;
        mTaskCount = maximum;
        mReporter = Executors.newSingleThreadScheduledExecutor(r -> new Thread(r, "ProgressTracker"));
        mReporter.scheduleAtFixedRate(this::reportProgress, 0, reportFrequencySeconds, TimeUnit.SECONDS);
        mStartTimeNanos = System.nanoTime();
    }

    public void increment()
    {
        increment(1);
    }

    public void increment(final int progress)
    {
        final int newValue = mTasksComplete.addAndGet(progress);
        if (newValue == mTaskCount)
        {
            reportProgress();
            mReporter.shutdownNow();
        }
    }

    public void reportProgress()
    {
        /*
        final int tasksCompleted = mTasksComplete.get();
        final double percentComplete = ((double) tasksCompleted / mTaskCount) * 100.0;

        final long currentTime = System.nanoTime();
        final long elapsedNanos = currentTime - mStartTimeNanos;

        final long nanosPerOperation = updateNanosPerOperation(currentTime, tasksCompleted);

        // Until we're near the end, artificially inflate the remaining time
        final int fudgeStopPercent = 66;
        final double maxFudgeFactor = 1; // Double the estimated time initially
        final double fudgeFactor = percentComplete > fudgeStopPercent
                ? 1
                : 1 + (((fudgeStopPercent - percentComplete) / fudgeStopPercent) * maxFudgeFactor);

        final long nanosRemaining = (long) ((mTaskCount - tasksCompleted) * (nanosPerOperation * fudgeFactor));

        LOGGER.info("{}{}% complete, {}/{} ({}/op) {} elapsed, ~{} remaining (this step)",
                mDescription == null ? "" : mDescription + ": ", String.format("%.2f", percentComplete), tasksCompleted, mTaskCount,
                StringUtils.formatNanos(nanosPerOperation), StringUtils.formatNanos(elapsedNanos), StringUtils.formatNanos(nanosRemaining));
        */
    }

    private long updateNanosPerOperation(final long currentTimeNanos, final int currentTasksCompleted)
    {
        final long previousTimeNanos;
        final int previousTasks;
        if (mPreviousTimeNanos.isEmpty())
        {
            previousTimeNanos = mStartTimeNanos;
            previousTasks = 0;
        }
        else
        {
            assert !mPreviousTasks.isEmpty();

            previousTimeNanos = mPreviousTimeNanos.peekFirst();
            previousTasks = mPreviousTasks.peekFirst();
        }

        mPreviousTimeNanos.addLast(currentTimeNanos);
        mPreviousTasks.addLast(currentTasksCompleted);
        if (mPreviousTimeNanos.size() > 5)
        {
            mPreviousTimeNanos.removeFirst();
            mPreviousTasks.removeFirst();
        }

        final long elapsedNanos = currentTimeNanos - previousTimeNanos;
        final int completedTasks = currentTasksCompleted - previousTasks;
        return completedTasks == 0 ? 0 : elapsedNanos / completedTasks;
    }

    @Override
    public void close()
    {
        mReporter.shutdownNow();
    }
}
