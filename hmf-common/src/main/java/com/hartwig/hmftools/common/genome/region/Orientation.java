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

    // byte methods for backwards compatibility
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

        throw new IllegalArgumentException("Invalid orientation: " + orientation);
    }

    public static Orientation fromByteStr(final String orientationStr) { return fromByte(Byte.parseByte(orientationStr)); }

    public boolean equalsByte(final byte orientation) { return (orientation == ORIENT_FWD) == (this == FORWARD); }
    public boolean isForward(final byte orientation) { return orientation == ORIENT_FWD; }
    public boolean isReverse(final byte orientation) { return orientation == ORIENT_FWD; }

    public static Orientation fromChar(char orientation)
    {
        switch(orientation)
        {
            case ORIENT_POS_CHAR:
                return Orientation.FORWARD;
            case ORIENT_NEG_CHAR:
                return Orientation.REVERSE;
        }

        throw new IllegalArgumentException("Invalid orientation: " + orientation);
    }

    public char asChar() { return this == FORWARD ? ORIENT_POS_CHAR : ORIENT_NEG_CHAR; }

    public Orientation opposite() { return opposite(this); }
    public static Orientation opposite(final Orientation orientation) { return orientation == FORWARD ? REVERSE : FORWARD; }

    public String toString() { return String.valueOf(asByte()); }
}
