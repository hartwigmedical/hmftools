package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import java.util.Set;
import java.util.StringJoiner;

import com.beust.jcommander.internal.Sets;

public enum ReadFilterType
{
    MIN_ALIGN_MATCH(1, 0, "Aligned base count"),
    MIN_MAP_QUAL(2, 1, "Min map quality"),
    INSERT_MAP_OVERLAP(4, 2,"Insert size vs aligned bases"),
    SOFT_CLIP_LENGTH(8, 3, "Soft-clip length"),
    SOFT_CLIP_BASE_QUAL(16, 4,  "Soft-clip base qual");
    // NO_SUPPORTING_READ(32, 5, "Read support");

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

    public String label()
    {
        return name().toLowerCase().replace('_', ' ');
    }

    public String description()
    {
        return mDescription;
    }

    public static ReadFilterType valueOf(int flag)
    {
        ReadFilterType[] allTypes = values();

        for(int type = 0; type < allTypes.length; ++type)
        {
            ReadFilterType f = allTypes[type];

            if(flag == f.mFlag)
                return f;
        }

        return null;
    }

    public static ReadFilterType findByName(String flag)
    {
        ReadFilterType[] allValues = values();

        for(int type = 0; type < allValues.length; ++type)
        {
            ReadFilterType f = allValues[type];
            if(f.name().equals(flag))
                return f;
        }

        return null;
    }

    public boolean isSet(int flag) { return (mFlag & flag) != 0; }

    public boolean isUnset(int flag) { return !isSet(flag); } // this.mFlags &= ~bit;

    public static Set<ReadFilterType> getFlags(int flag)
    {
        Set<ReadFilterType> set = Sets.newHashSet();
        ReadFilterType[] allValues = values();

        for(int type = 0; type < allValues.length; ++type)
        {
            ReadFilterType f = allValues[type];
            if(f.isSet(flag))
                set.add(f);
        }

        return set;
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
