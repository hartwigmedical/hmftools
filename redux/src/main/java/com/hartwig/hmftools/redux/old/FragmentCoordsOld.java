package com.hartwig.hmftools.redux.old;

import static java.lang.Math.abs;
import static java.lang.String.format;

public class FragmentCoordsOld
{
    public final String Key;
    public final int InitialPosition; // negative/reverse strand positions are negated
    public final boolean IsForward; // forward = F1R2, reverse is F2R1 - relates to collapsing and dual-strand classification
    public final boolean Incomplete;

    public static final FragmentCoordsOld NO_COORDS = new FragmentCoordsOld("", 0, true, true);

    public static final char FRAGMENT_REVERSED_ID = 'N';

    public FragmentCoordsOld(final String key, final int initialPosition, boolean isForward)
    {
        this(key, initialPosition, isForward, false);
    }

    public FragmentCoordsOld(final String key, int initialPosition, boolean isForward, boolean incomplete)
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

    public boolean matches(final FragmentCoordsOld other, boolean requireOrientation)
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
