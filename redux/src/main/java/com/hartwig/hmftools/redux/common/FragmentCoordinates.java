package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

public class FragmentCoordinates
{
    public final String Key;
    public final int InitialPosition; // negative/reverse strand positions are negated
    public final boolean IsForward; // forward = F1R2, reverse is F2R1 - relates to collapsing and dual-strand classification
    public final boolean Incomplete;

    public static final FragmentCoordinates NO_COORDS = new FragmentCoordinates("", 0, true, true);

    public static final char FRAGMENT_REVERSED_ID = 'N';

    public FragmentCoordinates(final String key, final int initialPosition, boolean isForward)
    {
        this(key, initialPosition, isForward, false);
    }

    public FragmentCoordinates(final String key, int initialPosition, boolean isForward, boolean incomplete)
    {
        Key = key;
        InitialPosition = initialPosition;
        IsForward = isForward;
        Incomplete = incomplete;
    }

    public String keyOriented() { return IsForward ? Key : format("%s_%c", Key, FRAGMENT_REVERSED_ID); }

    public String toString()
    {
        return Incomplete ? format("%s incomplete", keyOriented()) : keyOriented();
    }

    public boolean matches(final FragmentCoordinates other, boolean requireOrientation)
    {
        return Key.equals(other.Key) && (!requireOrientation || IsForward == other.IsForward);
    }

    public static String formCoordinate(final String chromosome, final int position, final boolean isForward)
    {
        return isForward ? format("%s_%d", chromosome, position) : format("%s_%d_R", chromosome, position);
    }

    public static String formKey(final String firstCoord, final String secondCoord)
    {
        return format("%s_%s", firstCoord, secondCoord);
    }

    public static String formKey(final String readCoord, final int insertSize)
    {
        return format("%s_%d", readCoord, abs(insertSize));
    }
}
