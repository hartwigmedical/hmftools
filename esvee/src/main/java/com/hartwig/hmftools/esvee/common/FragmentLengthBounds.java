package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

public class FragmentLengthBounds
{
    public final int UpperBound;
    public final int LowerBound;

    public static final FragmentLengthBounds INVALID = new FragmentLengthBounds(0, 0);

    public FragmentLengthBounds(final int lowerBound, final int upperBound)
    {
        UpperBound = upperBound;
        LowerBound = lowerBound;
    }

    public boolean isValid() { return LowerBound > 0 && UpperBound > LowerBound; }

    public String toString() { return format("lower(%d) upper(%d)", LowerBound, UpperBound); }
}
