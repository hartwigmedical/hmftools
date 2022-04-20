package com.hartwig.hmftools.isofox.adjusts;

public class FragmentSize
{
    public final int Length;
    public int Frequency;

    public FragmentSize(final int length, final int frequency)
    {
        Length = length;
        Frequency = frequency;
    }

    public String toString() { return String.format("length(%d) freq(%d)", Length, Frequency); }
}
