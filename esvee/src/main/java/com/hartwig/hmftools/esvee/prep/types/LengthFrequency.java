package com.hartwig.hmftools.esvee.prep.types;

public class LengthFrequency
{
    public final int Length;
    public int Frequency;

    public LengthFrequency(final int length, final int frequency)
    {
        Length = length;
        Frequency = frequency;
    }

    public String toString()
    {
        return String.format("length(%d) freq(%d)", Length, Frequency);
    }
}
