package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

public class FragmentCoordinates
{
    public final String Key;
    public final int InitialPosition; // negative/reverse strand positions are negated
    public final boolean Incomplete;

    public static final FragmentCoordinates NO_COORDS = new FragmentCoordinates("", 0, true);

    public FragmentCoordinates(final String key, final int initialPosition)
    {
        this(key, initialPosition, false);
    }

    public FragmentCoordinates(final String key, final int initialPosition, boolean incomplete)
    {
        Key = key;
        InitialPosition = initialPosition;
        Incomplete = incomplete;
    }

    public String toString()
    {
        return Incomplete ? format("%s incomplete", Key) : Key;
    }

    public static String formCoordinate(final String chromosome, final int position, final boolean isForward)
    {
        return isForward ? format("%s_%d", chromosome, position) : format("%s_%d_R", chromosome, position);
    }
}
