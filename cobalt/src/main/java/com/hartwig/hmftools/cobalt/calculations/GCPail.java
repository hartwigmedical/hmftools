package com.hartwig.hmftools.cobalt.calculations;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class GCPail
{
    public static int bucketIndex(final double gcContent)
    {
        return (int) Math.round(gcContent * 100);
    }

    public final int mGC;
    private final DescriptiveStatistics mStatistics = new DescriptiveStatistics();

    public GCPail(final int mGC)
    {
        Preconditions.checkArgument(mGC >= 0);
        Preconditions.checkArgument(mGC <= 100);
        this.mGC = mGC;
    }

    public double median()
    {
        if(mStatistics.getN() == 0)
        {
            return 0;
        }
        return mStatistics.getPercentile(50);
    }

    public void addReading(final double value)
    {
        mStatistics.addValue(value);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final GCPail gcPail = (GCPail) o;
        return mGC == gcPail.mGC;
    }

    @Override
    public int hashCode()
    {
        return Objects.hashCode(mGC);
    }

    @Override
    public String toString()
    {
        return "GCPail{" +
                "mGC=" + mGC +
                '}';
    }
}
