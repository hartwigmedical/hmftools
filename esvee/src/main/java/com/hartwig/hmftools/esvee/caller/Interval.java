package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.abs;

import java.util.List;

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

    public static Interval fromCiposTag(List<Integer> ciposValues)
    {
        if(ciposValues == null || ciposValues.size() != 2)
            return new Interval();

        return new Interval(ciposValues.get(0), ciposValues.get(1));
    }

    public int length() { return abs(End - Start); }

    public boolean matches(final Interval other)
    {
        return Start == other.Start && End == other.End;
    }

    public String toString() { return String.format("[%d,%d]", Start, End); }
}
