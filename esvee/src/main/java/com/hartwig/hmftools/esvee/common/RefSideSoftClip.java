package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class RefSideSoftClip
{
    public final int Position;
    public final byte Orientation;

    private final List<Read> mReads;
    private int mMaxLength;

    public RefSideSoftClip(final int position, final byte orientation, final Read read, final int readSoftClipLength)
    {
        Position = position;
        Orientation = orientation;
        mReads = Lists.newArrayList(read);
        mMaxLength = readSoftClipLength;
    }

    public List<Read> reads() { return mReads; }
    public int maxLength() { return mMaxLength; }
    public int readCount() { return mReads.size(); }

    public void addRead(final Read read, final int readSoftClipLength)
    {
        mReads.add(read);
        mMaxLength = max(mMaxLength, readSoftClipLength);
    }

    public String toString() { return format("%d len(%d) reads(%d)", Position, mMaxLength, mReads.size()); }
}
