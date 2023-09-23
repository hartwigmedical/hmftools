package com.hartwig.hmftools.sieve.annotate;

import org.jetbrains.annotations.NotNull;

public class RepeatMasker
{
    public static final String CSV_HEADER = "RepeatType,RepeatInfo,RepeatPosStart,RepeatPosEnd,Count,OtherInfo";
    public static final String EMPTY_CSV_FRAGMENT = "NA,NA,NA,NA,NA,NA";

    private final String mRepeatType;
    private final String mRepeatInfo;
    private final int mRepeatPosStart;
    private final int mRepeatPosEnd;
    private final int mCount;
    private final String mOtherInfo;

    public RepeatMasker(@NotNull final String repeatType, @NotNull final String repeatInfo, @NotNull final int repeatPosStart,
            @NotNull final int repeatPosEnd, @NotNull final int count, @NotNull final String otherInfo)
    {
        mRepeatType = repeatType;
        mRepeatInfo = repeatInfo;
        mRepeatPosStart = repeatPosStart;
        mRepeatPosEnd = repeatPosEnd;
        mCount = count;
        mOtherInfo = otherInfo;
    }

    public String getCSVFragment()
    {
        return mRepeatType
                + ','
                + mRepeatInfo
                + ','
                + mRepeatPosStart
                + ','
                + mRepeatPosEnd
                + ','
                + mCount
                + ','
                + mOtherInfo;
    }

    public String getRepeatType()
    {
        return mRepeatType;
    }

    public String getRepeatInfo()
    {
        return mRepeatInfo;
    }

    public int getRepeatPosStart()
    {
        return mRepeatPosStart;
    }

    public int getRepeatPosEnd()
    {
        return mRepeatPosEnd;
    }

    public int getCount()
    {
        return mCount;
    }

    public String getOtherInfo()
    {
        return mOtherInfo;
    }

    @Override
    public String toString()
    {
        return "RepeatMasker{" +
                "mepeatType='" + mRepeatType + '\'' +
                ", RepeatInfo='" + mRepeatInfo + '\'' +
                ", RepeatPosStart=" + mRepeatPosStart +
                ", RepeatPosEnd=" + mRepeatPosEnd +
                ", Count=" + mCount +
                ", OtherInfo='" + mOtherInfo + '\'' +
                '}';
    }
}
