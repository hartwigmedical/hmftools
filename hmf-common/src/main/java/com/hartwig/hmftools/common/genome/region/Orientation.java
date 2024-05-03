package com.hartwig.hmftools.common.genome.region;

public enum Orientation
{
    FORWARD,
    REVERSE;

    public static final byte ORIENT_FWD = 1;
    public static final byte ORIENT_REV = -1;

    public static final char ORIENT_POS_CHAR = '+';
    public static final char ORIENT_NEG_CHAR = '-';

    public boolean isForward() { return this == FORWARD; }
    public boolean isReverse() { return this == REVERSE; }

    public byte asByte() { return this == FORWARD ? ORIENT_FWD : ORIENT_REV; }

    public static Orientation fromByte(byte orientation)
    {
        switch(orientation)
        {
            case 1:
                return Orientation.FORWARD;
            case -1:
                return Orientation.REVERSE;
        }

        throw new IllegalArgumentException("Invalid direction: " + orientation);
    }

    public static Orientation fromStr(final String orientation) { return fromChar(orientation.charAt(0)); }

    public static Orientation fromChar(char orientation)
    {
        switch(orientation)
        {
            case ORIENT_POS_CHAR:
                return Orientation.FORWARD;
            case ORIENT_NEG_CHAR:
                return Orientation.REVERSE;
        }

        throw new IllegalArgumentException("Invalid strand: " + orientation);
    }

    public char asChar() { return this == FORWARD ? ORIENT_POS_CHAR : ORIENT_NEG_CHAR; }

    public Orientation opposite() { return opposite(this); }
    public static Orientation opposite(final Orientation orientation) { return orientation == FORWARD ? REVERSE : FORWARD; }

    public String toString() { return String.valueOf(asByte()); }
}
