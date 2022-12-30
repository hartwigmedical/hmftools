package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

public class FragmentCoordinates
{
    public final String Key;
    public final int InitialPosition; // negative/reverse strand positions are negated

    public FragmentCoordinates(final String key, final int initialPosition)
    {
        Key = key;
        InitialPosition = initialPosition;
    }

    public String toString() { return Key; }

    public static String formCoordinate(final String chromosome, final int position, final boolean isForward)
    {
        return isForward ? format("%s_%d", chromosome, position) : format("%s_%d_R", chromosome, position);
    }
}
