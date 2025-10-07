package com.hartwig.hmftools.redux.common;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.redux.duplicate.FragmentCoords;

import htsjdk.samtools.SAMRecord;

public class ReadInfo
{
    private final SAMRecord mRead;
    private final FragmentCoords mCoordinates;

    public ReadInfo(final SAMRecord read, final FragmentCoords fragCoords)
    {
        mRead = read;
        mCoordinates = fragCoords;
    }

    public String id() { return mRead.getReadName(); }
    public SAMRecord read() { return mRead; }

    public FragmentCoords coordinates() { return mCoordinates; }

    public String toString()
    {
        return String.format("read(%s:%d %s) coords(%s)",
                mRead.getReferenceName(), mRead.getAlignmentStart(), id(), mCoordinates);
    }

    public static String readToString(final SAMRecord read) { return SamRecordUtils.readToString(read); }
}
