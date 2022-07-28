package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

import java.util.StringJoiner;

public enum ReadFilterType
{
    MIN_ALIGN_MATCH(1, 0, "Aligned base count"),
    MIN_MAP_QUAL(2, 1, "Min map quality"),
    INSERT_MAP_OVERLAP(4, 2,"Insert size vs aligned bases"),
    SOFT_CLIP_LENGTH(8, 3, "Soft-clip length"),
    SOFT_CLIP_BASE_QUAL(16, 4,  "Soft-clip insufficient high base qual"),
    BREAK_IN_REPEAT(32, 5,  "Repeat break"),
    POLY_G_SC(64, 6,  "Poly-G"),
    SOFT_CLIP_LOW_BASE_QUAL(128, 6,  "Soft-clip excessive low base qual");

    private final int mFlag;
    private final int mIndex;
    private final String mDescription;

    ReadFilterType(int flag, int index, String description)
    {
        mFlag = flag;
        mIndex = index;
        mDescription = description;
    }

    public int flag()
    {
        return mFlag;
    }
    public int index() { return mIndex; }
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
}
