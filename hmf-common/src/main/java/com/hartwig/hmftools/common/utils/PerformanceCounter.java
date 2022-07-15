package com.hartwig.hmftools.common.utils;

import static java.lang.Math.floor;
import static java.lang.Math.max;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PerformanceCounter
{
    private static final Logger LOGGER = LogManager.getLogger(PerformanceCounter.class);

    private final String mName;
    // time values are in seconds
    private double mTotalTime;
    private double mMaxTime;
    private boolean mIsRunning;
    private boolean mIsPaused;

    private long mStartTime;
    private long mPausedTime; // accumulates interval times when timer is paused
    private final List<Double> mTimes;
    private final List<String> mTimeNames;
    private String mCurrentIntervalName;
    private boolean mSortTimes;

    private static final double NANOS_IN_SECOND = 1000000000;

    public PerformanceCounter(final String name)
    {
        mName = name;
        mIsRunning = false;
        mStartTime = 0;
        mPausedTime = 0;
        mTotalTime = 0;
        mMaxTime = 0;
        mTimes = Lists.newArrayList();
        mTimeNames = Lists.newArrayList();
        mSortTimes = false;
        mCurrentIntervalName = null;
    }

    public void reset()
    {
        mIsRunning = false;
        mStartTime = 0;
        mPausedTime = 0;
        mTotalTime = 0;
        mMaxTime = 0;
        mTimes.clear();
        mTimeNames.clear();
        mCurrentIntervalName = null;
    }

    public void setSortTimes(boolean toggle) { mSortTimes = toggle; }

    public final String getName() {
        return mName;
    }

    public void start()
    {
        start(null);
    }

    public void start(final String intervalName)
    {
        mIsRunning = true;
        mIsPaused = false;
        mPausedTime = 0;
        mStartTime = System.nanoTime();

        mCurrentIntervalName = intervalName;
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

        mMaxTime = max(sampleTimeSeconds, mMaxTime);
        mTotalTime += sampleTimeSeconds;

        if(!mSortTimes)
        {
            mTimes.add(sampleTimeSeconds);

        }
        else
        {
            addTimeInOrder(sampleTimeSeconds, mCurrentIntervalName);
        }
    }

    private void addTimeInOrder(double time, final String name)
    {
        int index = 0;
        while(index < mTimes.size())
        {
            if(time < mTimes.get(index))
                break;
            else
                ++index;
        }

        mTimes.add(index, time);

        if(name != null)
            mTimeNames.add(index, name);
    }

    public boolean isRunning() {
        return mIsRunning;
    }

    public int getSampleCount() {
        return mTimes.size();
    }

    public final List<Double> getTimes() {
        return mTimes;
    }
    public final List<String> getTimeNames() { return mTimeNames; }

    // time values are in seconds
    public double getTotalTime() {
        return mTotalTime;
    }

    public double getMaxTime() {
        return mMaxTime;
    }

    public double getAvgTime() {
        return mTimes.isEmpty() ? 0 : mTotalTime / mTimes.size();
    }

    public double getMedianTime()
    {
        if(!mSortTimes || mTimes.isEmpty())
            return 0;

        int medianIndex = (int)floor(mTimes.size()/2D);
        return mTimes.get(medianIndex);
    }

    public double getLastTime() { return !mTimes.isEmpty() ? mTimes.get(mTimes.size() - 1) : 0; }

    public void logStats()
    {
        logStats(false);
    }

    public void logStats(boolean logIntervals)
    {
        if(mTimes.isEmpty())
            return;

        if(mTimes.size() > 1)
        {
            String avgMed =  String.format("avg=%.3f", getAvgTime());

            if(mSortTimes)
                avgMed += String.format(" med=%.3f", getMedianTime());

            LOGGER.info(String.format("PerfStats: name(%s) intervals(%d) total(%.3f %s max=%.3f)",
                    mName, getSampleCount(), getTotalTime(), avgMed, getMaxTime()));
        }
        else
        {
            LOGGER.info(String.format("PerfStats: name(%s) intervals(%d) total(%.3f)",
                    mName, getSampleCount(), getTotalTime()));
        }

        if(logIntervals && mTimes.size() > 1)
        {
            // log the individual interval data
            for(int i = 0; i < mTimes.size(); ++i)
            {
                if(mTimes.size() == mTimeNames.size())
                {
                    LOGGER.info(String.format("PerfStats: interval(%s) time(%.3f)", mTimeNames.get(i), mTimes.get(i)));
                }
                else
                {
                    LOGGER.info(String.format("PerfStats: interval %d: time(%.3f)", i, mTimes.get(i)));
                }
            }
        }
    }

    public void merge(final PerformanceCounter other)
    {
        if(!mSortTimes)
        {
            mTimes.addAll(other.getTimes());
            mTimeNames.addAll(other.getTimeNames());
        }
        else
        {
            final List<String> otherTimeNames = other.getTimeNames();
            for(int i = 0; i < other.getTimes().size(); ++i)
            {
                addTimeInOrder(other.getTimes().get(i), otherTimeNames.isEmpty() ? null : otherTimeNames.get(i));
            }
        }

        mTotalTime += other.getTotalTime();
        mMaxTime = max(mMaxTime, other.getMaxTime());
    }

    public static double nanosToSeconds(long nanoStartTime, long nanosEndTime)
    {
        return (nanosEndTime - nanoStartTime) / NANOS_IN_SECOND;
    }
}
