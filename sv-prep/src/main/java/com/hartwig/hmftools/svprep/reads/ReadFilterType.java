package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

import java.util.EnumSet;
import java.util.StringJoiner;

public enum ReadFilterType
{
    /** Too few bases matched */
    MIN_ALIGN_MATCH(1, "Aligned base count"),

    /** Map quality was too low */
    MIN_MAP_QUAL(2, "Min map quality"),

    /**  */
    INSERT_MAP_OVERLAP(4, "Insert size vs aligned bases"),

    /** Soft clip was too short */
    SOFT_CLIP_LENGTH(8, "Soft-clip length"),

    /** Soft clipped area average base qual was too low */
    SOFT_CLIP_BASE_QUAL(16,  "Soft-clip insufficient high base qual"),

    /** Many of the same bases have been insterted */
    BREAK_IN_REPEAT(32, "Repeat break"),

    /** Too many Gs or Cs were present */
    POLY_G_SC(64, "Poly-G"),

    /** Soft clipped area had too many bases of low quality */
    SOFT_CLIP_LOW_BASE_QUAL(128, "Soft-clip excessive low base qual");

    private final int mFlag;
    private final String mDescription;

    ReadFilterType(int flag, String description)
    {
        mFlag = flag;
        mDescription = description;
    }

    public int flag()
    {
        return mFlag;
    }
    public int index() { return ordinal(); }
    public String description()
    {
        return mDescription;
    }

    public boolean isSet(int flag) { return (mFlag & flag) != 0; }
    public static boolean isSet(int flags, ReadFilterType flag) { return (flags & flag.flag()) != 0; }
    public boolean isUnset(int flag) { return !isSet(flag); } // this.mFlags &= ~bit;

    public static int set(int flags, ReadFilterType flag)
    {
        flags |= flag.flag();
        return flags;
    }

    public static int unset(int flags, ReadFilterType flag)
    {
        flags &= ~flag.flag();
        return flags;
    }

    public static String filterCountsToString(final int[] filterCounts)
    {
        StringJoiner sj = new StringJoiner(", ");

        for(ReadFilterType type : ReadFilterType.values())
        {
            sj.add(format("%s=%d", type, filterCounts[type.index()]));
        }

        return sj.toString();
    }

    public static int toFlags(EnumSet<ReadFilterType> readFilterTypes) {
        int flags = 0;
        for (final ReadFilterType type : readFilterTypes) {
            flags |= type.flag();
        }
        return flags;
    }
}
