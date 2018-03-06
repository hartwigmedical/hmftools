package com.hartwig.hmftools.svannotation.analysis;


public class SvPON {

    final private String mChrStart;
    final private String mChrEnd;
    final private long mPosStart;
    final private long mPosEnd;
    final private byte mOrientStart;
    final private byte mOrientEnd;
    final private String mType;
    final private int mCount;

    public SvPON(
        final String chrStart,
        final String chrEnd,
        final long posStart,
        final long posEnd,
        final byte orientStart,
        final byte orientEnd,
        final String type,
        final int count)
    {
        mChrStart = chrStart;
        mChrEnd = chrEnd;
        mPosStart = posStart;
        mPosEnd = posEnd;
        mOrientStart = orientStart;
        mOrientEnd = orientEnd;
        mType = type;
        mCount = count;
    }

    public String chrStart() { return mChrStart; }
    public String chrEnd() { return mChrEnd; }
    public long posStart() { return mPosStart; }
    public long posEnd() { return mPosEnd; }
    public byte orientStart() { return mOrientStart; }
    public byte orientEnd() { return mOrientEnd; }
    public String type() { return mType; }
    public int count() { return mCount; }
}
