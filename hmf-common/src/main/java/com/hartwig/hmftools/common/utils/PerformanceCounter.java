package com.hartwig.hmftools.common.utils;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PerformanceCounter
{
    private static final Logger LOGGER = LogManager.getLogger(PerformanceCounter.class);

    private final String mName;
    private boolean mTrackTimes; // keeps each recorded time, for top-N logging and median calcs

    // time values are in seconds
    private double mTotalTime;
    private double mMaxTime;
    private boolean mIsRunning;
    private boolean mIsPaused;

    private long mStartTime;
    private long mPausedTime; // accumulates interval times when timer is paused
    private int mIntervalCount;
    private List<NamedTime> mNamedTimes;

    private double mLastTime;
    private String mCurrentIntervalName;

    private static final double NANOS_IN_SECOND = 1000000000;

    public PerformanceCounter(final String name)
    {
        mName = name;
        mTrackTimes = true;
        mIsRunning = false;
        mStartTime = 0;
        mPausedTime = 0;
        mTotalTime = 0;
        mMaxTime = 0;
        mNamedTimes = null;
        mIntervalCount = 0;
        mLastTime = 0;
        mCurrentIntervalName = null;
    }

    public void reset()
    {
        mIsRunning = false;
        mStartTime = 0;
        mPausedTime = 0;
        mTotalTime = 0;
        mMaxTime = 0;
        mIntervalCount = 0;
        mLastTime = 0;
        mCurrentIntervalName = null;

        if(mNamedTimes != null)
            mNamedTimes.clear();
    }

    public final String getName() {
        return mName;
    }
    public List<NamedTime> getNamedTimes() { return mNamedTimes; }

    public void setTrackTimes() { mTrackTimes = true; }

    public void start()
    {
        if(mTrackTimes)
            start(String.valueOf(mIntervalCount));
        else
            start(null);
    }

    public void start(final String intervalName)
    {
        mIsRunning = true;
        mIsPaused = false;
        mPausedTime = 0;
        mStartTime = System.nanoTime();

        if(intervalName != null)
        {
            mCurrentIntervalName = intervalName;

            if(mNamedTimes == null)
                mNamedTimes = Lists.newArrayList();
        }
    }

    public void pause()
    {
        if(!mIsRunning)
            return;

        mIsPaused = true;
        mPausedTime += System.nanoTime() - mStartTime;
    }

    public void resume()
    {
        if(!mIsPaused)
            return;

        mIsPaused = false;
        mStartTime = System.nanoTime();
    }

    public void stop()
    {
        if (!mIsRunning)
            return;

        mIsRunning = false;
        mIsPaused = false;

        long sampleTime = System.nanoTime() - mStartTime;
        sampleTime += mPausedTime;
        double sampleTimeSeconds = sampleTime / NANOS_IN_SECOND;

        mLastTime = sampleTimeSeconds;
        mMaxTime = max(sampleTimeSeconds, mMaxTime);
        mTotalTime += sampleTimeSeconds;
        ++mIntervalCount;

        if(mCurrentIntervalName != null)
            mNamedTimes.add(new NamedTime(sampleTimeSeconds, mCurrentIntervalName));
    }

    public boolean isRunning() { return mIsRunning; }

    public int getIntervalCount() { return mIntervalCount; }

    // time values are in seconds
    public double getTotalTime() { return mTotalTime; }
    public double getMaxTime() { return mMaxTime; }
    public double getLastTime() { return mLastTime; }

    public double getAvgTime() { return mIntervalCount > 0 ? mTotalTime / (double)mIntervalCount : 0; }

    public double getMedianTime()
    {
        if(mNamedTimes == null || mNamedTimes.isEmpty())
            return 0;

        Collections.sort(mNamedTimes, new TimeComparator());

        int medianIndex = (int)floor(mNamedTimes.size()/2D);
        return mNamedTimes.get(medianIndex).Time;
    }

    public void logStats()
    {
        if(mIntervalCount == 0)
            return;

        LOGGER.info(String.format("PerfStat(%s) intervals(%d) total(%.3f avg=%.3f max=%.3f)",
                mName, getIntervalCount(), getTotalTime(), getAvgTime(), getMaxTime()));
    }

    public void logIntervalStats() { logIntervalStats(10); }

    public void logIntervalStats(int topN)
    {
        if(mIntervalCount == 0)
            return;

        if(mNamedTimes == null)
        {
            logStats();
            return;
        }

        LOGGER.info(String.format("PerfStat(%s) intervals(%d) total(%.3f avg=%.3f med=%.3f max=%.3f)",
                mName, getIntervalCount(), getTotalTime(), getAvgTime(), getMedianTime(), getMaxTime()));

        // median call above will have sorted the times

        // log the individual interval data
        int maxTimes = topN > 0 ? min(mNamedTimes.size(), topN) : mNamedTimes.size();
        for(int i = 0; i < maxTimes; ++i)
        {
            LOGGER.info(String.format("PerfStats(%s) interval(%s) time(%.3f)", mName, mNamedTimes.get(i).Name, mNamedTimes.get(i).Time));
        }
    }

    public void merge(final PerformanceCounter other)
    {
        // assumes both counters are stopped / complete
        mTotalTime += other.getTotalTime();
        mMaxTime = max(mMaxTime, other.getMaxTime());
        mIntervalCount += other.getIntervalCount();

        if(other.getNamedTimes() != null)
        {
            if(mNamedTimes == null)
                mNamedTimes = Lists.newArrayList();

            mNamedTimes.addAll(other.getNamedTimes());
        }
    }

    private class NamedTime
    {
        public final double Time;
        public final String Name;

        public NamedTime(final double time, final String name)
        {
            Time = time;
            Name = name;
        }
    }

    public static class TimeComparator implements Comparator<NamedTime>
    {
        public int compare(final NamedTime first, final NamedTime second)
        {
            return first.Time < second.Time ? 1 : -1;
        }
    }

    public static double nanosToSeconds(long nanoStartTime, long nanosEndTime)
    {
        return (nanosEndTime - nanoStartTime) / NANOS_IN_SECOND;
    }
}
