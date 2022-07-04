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
    public final List<ReadGroup> JunctionGroups; // with a read matching the junction
    public final List<ReadRecord> SupportingReads;
    public final List<RemoteJunction> RemoteJunctions;

    private boolean mHotspot;

    public JunctionData(final int position, final byte orientation, final ReadRecord read)
    {
        Position = position;
        Orientation = orientation;

        JunctionGroups = Lists.newArrayList();
        SupportingReads = Lists.newArrayList();
        RemoteJunctions = Lists.newArrayList();
        InitialRead = read;

        mHotspot = false;
    }

    public int exactFragmentCount() { return JunctionGroups.size(); }
    public int supportingReadCount() { return SupportingReads.size(); }
    public int totalSupport() { return exactFragmentCount() + supportingReadCount(); }
    public boolean hotspot() { return mHotspot; }
    public void markHotspot() { mHotspot = true; }

    public String toString()
    {
        return format("loc(%d:%d) frags(exact=%d supp=%d)", Position, Orientation, exactFragmentCount(), supportingReadCount());
    }

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
