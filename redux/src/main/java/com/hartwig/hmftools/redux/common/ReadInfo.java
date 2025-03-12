package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import htsjdk.samtools.SAMRecord;

public class ReadInfo
{
    private final SAMRecord mRead;

    // with duplicate group collapsing mCoordinates may not match mRead, see mPreCollapsedCoordinates
    private FragmentCoords mCoordinates;
    private FragmentCoords mPreCollapsedCoordinates;

    public ReadInfo(final SAMRecord read, final FragmentCoords fragCoords)
    {
        mRead = read;
        mCoordinates = fragCoords;
        mPreCollapsedCoordinates = null;
    }

    public String id() { return mRead.getReadName(); }
    public SAMRecord read() { return mRead; }

    public FragmentCoords coordinates() { return mCoordinates; }
    public FragmentCoords preCollapsedCoordinates() { return mPreCollapsedCoordinates == null ? mCoordinates : mPreCollapsedCoordinates; }

    public void updateCoordinates(final FragmentCoords coords)
    {
        if(mCoordinates.equals(coords))
            return;

        if(mPreCollapsedCoordinates == null)
            mPreCollapsedCoordinates = mCoordinates;

        mCoordinates = coords;
    }

    public String toString()
    {
        return String.format("read(%s:%d %s) coords(%s)",
                mRead.getReferenceName(), mRead.getAlignmentStart(), id(), mCoordinates);
    }

    public static String readToString(final SAMRecord read)
    {
        if(read.getReadPairedFlag())
        {
            return format("id(%s) coords(%s:%d-%d) isPaired(true) cigar(%s) mate(%s:%d) flags(%d)",
                    read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                    read.getCigarString(), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getFlags());
        }

        return format("id(%s) coords(%s:%d-%d) isPaired(false) cigar(%s) flags(%d)",
                read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                read.getCigarString(), read.getFlags());
    }
}
