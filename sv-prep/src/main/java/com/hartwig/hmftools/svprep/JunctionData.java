package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;

public class JunctionData
{
    public final int Position;
    public final byte Orientation;

    public final ReadRecord InitialRead;
    //public final List<ReadRecord> Reads; // matching the junction
    public final List<RemoteJunction> RemoteJunctions;

    public int ExactFragments;
    public int SupportReads;
    private boolean mHotspot;

    public JunctionData(final int position, final byte orientation, final ReadRecord read)
    {
        Position = position;
        Orientation = orientation;

        // Reads = Lists.newArrayList();
        RemoteJunctions = Lists.newArrayList();
        InitialRead = read;

        ExactFragments = 0;
        SupportReads = 0;
        mHotspot = false;
    }

    public int totalSupport() { return ExactFragments + SupportReads; }
    public boolean hotspot() { return mHotspot; }
    public void markHotspot() { mHotspot = true; }

    public String toString() { return format("loc(%d:%d) frags(exact=%d supp=%d)", Position, Orientation, ExactFragments, SupportReads); }

    public static class JunctionDataSorter implements Comparator<JunctionData>
    {
        // sorts by support descending
        public int compare(final JunctionData first, final JunctionData second)
        {
            if(first.totalSupport() != second.totalSupport())
                return first.totalSupport() < second.totalSupport() ? 1 : -1;

            return 0;
        }
    }

}
