package com.hartwig.hmftools.svprep.reads;

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
    public final List<ReadGroup> SupportingGroups;
    public final List<RemoteJunction> RemoteJunctions;

    private boolean mHotspot;

    public JunctionData(final int position, final byte orientation, final ReadRecord read)
    {
        Position = position;
        Orientation = orientation;

        JunctionGroups = Lists.newArrayList();
        SupportingGroups = Lists.newArrayList();
        RemoteJunctions = Lists.newArrayList();
        InitialRead = read;

        mHotspot = false;
    }

    public int exactFragmentCount() { return JunctionGroups.size(); }
    public int supportingReadCount() { return SupportingGroups.size(); }
    public int totalSupport() { return exactFragmentCount() + supportingReadCount(); }
    public boolean hotspot() { return mHotspot; }
    public void markHotspot() { mHotspot = true; }

    public String toString()
    {
        return format("loc(%d:%d) frags(exact=%d supp=%d) remotes(%d)",
                Position, Orientation, exactFragmentCount(), supportingReadCount(), RemoteJunctions.size());
    }
}
