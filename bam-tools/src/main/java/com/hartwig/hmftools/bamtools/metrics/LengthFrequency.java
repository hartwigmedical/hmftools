package com.hartwig.hmftools.bamtools.metrics;

public class LengthFrequency
{
    public final int Length;
    public long Frequency;

    public LengthFrequency(final int length, final long frequency)
    {
        Length = length;
        Frequency = frequency;
    }

    public String toString()
    {
        return String.format("length(%d) freq(%d)", Length, Frequency);
    }
}
