package com.hartwig.hmftools.redux.common;

import static com.hartwig.hmftools.redux.common.FragmentStatus.UNSET;

import htsjdk.samtools.SAMRecord;

public class ReadInfo
{
    private final SAMRecord mRead;
    private FragmentCoords mCoordinates;

    private boolean mReadsWritten;

    // duplicate read info
    private FragmentStatus mStatus;
    private String mUmi;

    public ReadInfo(final SAMRecord read, final FragmentCoords fragCoords)
    {
        mRead = read;
        mStatus = UNSET;
        mCoordinates = fragCoords;
        mReadsWritten = false;
        mUmi = null;
    }

    public String id() { return mRead.getReadName(); }
    public SAMRecord read() { return mRead; }

    public FragmentStatus status() { return mStatus; }
    public void setStatus(final FragmentStatus status) { mStatus = status; }

    public FragmentCoords coordinates() { return mCoordinates; }

    public String umi() { return mUmi; }
    public void setUmi(final String umi) { mUmi = umi; }

    public boolean readsWritten() { return mReadsWritten; }
    public void setReadWritten() { mReadsWritten = true; }

    public String toString()
    {
        return String.format("reads(%s:%d %s) status(%s) coords(%s)",
                mRead.getReferenceName(), mRead.getAlignmentStart(), id(), mStatus, mCoordinates);
    }
}
