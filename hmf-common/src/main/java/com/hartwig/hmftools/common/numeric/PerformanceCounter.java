package com.hartwig.hmftools.common.numeric;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PerformanceCounter {

    private final String mName;
    // time values are in seconds
    private double mTotalTime;
    private double mMaxTime;
    private boolean mIsRunning;

    private long mStartTime;
    private List<Double> mTimes;

    private static double NANOS_IN_SECOND = 1000000000;

    private static final Logger LOGGER = LogManager.getLogger(PerformanceCounter.class);

    public PerformanceCounter(final String name) {

        mName = name;
        mIsRunning = false;
        mStartTime = 0;
        mTotalTime = 0;
        mMaxTime = 0;
        mTimes = Lists.newArrayList();
    }

    public final String getName() { return mName; }

    public void start()
    {
        mIsRunning = true;
        mStartTime = System.nanoTime();
    }

    public void stop()
    {
        if(!mIsRunning)
            return;

        mIsRunning = false;

        long sampleTime = System.nanoTime() - mStartTime;
        double sampleTimeSeconds = sampleTime / NANOS_IN_SECOND;

        mMaxTime = Math.max(sampleTimeSeconds, mMaxTime);
        mTotalTime += sampleTimeSeconds;

        mTimes.add(sampleTimeSeconds);
    }

    public boolean isRunning() { return mIsRunning; }

    public int getSampleCount() { return mTimes.size(); }
    public final List<Double> getSamples() { return mTimes; }

    // time values are in seconds
    public double getTotalTime() { return mTotalTime; }
    public double getMaxTime() { return mMaxTime; }
    public double getAvgTime() { return mTimes.isEmpty() ? 0 : mTotalTime / mTimes.size(); }

    public void logStats()
    {
        LOGGER.info(String.format("name(%s) samples(%d) total(%.3f) avg(%.3f) max(%.3f)",
                mName, getSampleCount(), getTotalTime(), getAvgTime(), getMaxTime()));

    }
}
