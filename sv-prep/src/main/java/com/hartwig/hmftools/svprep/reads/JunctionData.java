package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

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

    private boolean mInternalIndel;
    private boolean mHotspot;
    private int mDepth;

    public JunctionData(final int position, final byte orientation, final ReadRecord read)
    {
        Position = position;
        Orientation = orientation;

        JunctionGroups = Lists.newArrayList();
        SupportingGroups = Lists.newArrayList();
        RemoteJunctions = Lists.newArrayList();
        InitialRead = read;

        mHotspot = false;
        mInternalIndel = false;
        mDepth = 0;
    }

    public boolean isExisting() { return InitialRead == null; }

    public int junctionFragmentCount() { return JunctionGroups.size(); }
    public int supportingFragmentCount() { return SupportingGroups.size(); }

    public boolean hotspot() { return mHotspot; }
    public void markHotspot() { mHotspot = true; }

    public boolean internalIndel() { return mInternalIndel; }
    public void markInternalIndel() { mInternalIndel = true; }

    public void setDepth(int depth) { mDepth = depth; }
    public int depth() { return mDepth; }

    public void addRemoteJunction(final RemoteJunction remoteJunction)
    {
        RemoteJunction matched = RemoteJunctions.stream().filter(x -> x.matches(remoteJunction)).findFirst().orElse(null);
        if(matched != null)
        {
            ++matched.Fragments;
            return;
        }

        RemoteJunctions.add(remoteJunction);
    }

    public String toString()
    {
        return format("loc(%d:%d) frags(junc=%d supp=%d) remotes(%d)",
                Position, Orientation, junctionFragmentCount(), supportingFragmentCount(), RemoteJunctions.size());
    }
}
