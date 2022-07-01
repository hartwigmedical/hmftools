package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import java.util.Comparator;

public class JunctionData
{
    public final int Position;
    public final byte Orientation;

    public final ReadRecord InitialRead;
    //public final List<ReadRecord> Reads; // matching the junction

    public int ExactReads;
    public int SupportReads;

    public JunctionData(final int position, final byte orientation, final ReadRecord read)
    {
        Position = position;
        Orientation = orientation;

        // Reads = Lists.newArrayList();
        InitialRead = read;

        ExactReads = 1;
        SupportReads = 0;
    }

    public int totalSupport() { return ExactReads + SupportReads; }

    public String toString() { return format("loc(%d:%d) reads(exact=%d supp=%d)", Position, Orientation, ExactReads, SupportReads); }

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
