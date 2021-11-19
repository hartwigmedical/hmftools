package com.hartwig.hmftools.gripss;

public class Interval
{
    public final int Start;
    public final int End;

    public Interval()
    {
        Start = 0;
        End = 0;
    }

    public Interval(final int start, final int end)
    {
        Start = start;
        End = end;
    }

    public boolean matches(final Interval other)
    {
        return Start == other.Start && End == other.End;
    }

    public String toString() { return String.format("[%d,%d]", Start, End); }
}
