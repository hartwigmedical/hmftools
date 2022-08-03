package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class JunctionData
{
    public final int Position;
    public final byte Orientation;

    public final ReadRecord InitialRead;
    public final List<ReadGroup> JunctionGroups; // with a read matching the junction
    public final List<ReadGroup> SupportingGroups;
    public final List<ReadGroup> ExactSupportGroups;
    public final List<RemoteJunction> RemoteJunctions;

    public final Map<ReadType,List<ReadRecord>> ReadTypeReads;

    private boolean mInternalIndel;
    private boolean mDiscordantGroup;
    private boolean mHotspot;
    private int mDepth;

    public JunctionData(final int position, final byte orientation, final ReadRecord read)
    {
        Position = position;
        Orientation = orientation;

        JunctionGroups = Lists.newArrayList();
        SupportingGroups = Lists.newArrayList();
        ExactSupportGroups = Lists.newArrayList();
        RemoteJunctions = Lists.newArrayList();
        ReadTypeReads = Maps.newHashMap();

        ReadTypeReads.put(ReadType.JUNCTION, Lists.newArrayList());
        ReadTypeReads.put(ReadType.SUPPORT, Lists.newArrayList());
        ReadTypeReads.put(ReadType.EXACT_SUPPORT, Lists.newArrayList());

        InitialRead = read;

        mHotspot = false;
        mInternalIndel = false;
        mDiscordantGroup = false;
        mDepth = 0;
    }

    public boolean isExisting() { return InitialRead == null; }

    public int junctionFragmentCount() { return JunctionGroups.size(); }
    public int supportingFragmentCount() { return SupportingGroups.size(); }
    public int exactSupportFragmentCount() { return ExactSupportGroups.size(); }
    public int totalFragmentCount() { return JunctionGroups.size() + SupportingGroups.size() + ExactSupportGroups.size(); }

    public boolean hotspot() { return mHotspot; }
    public void markHotspot() { mHotspot = true; }

    public boolean internalIndel() { return mInternalIndel; }
    public void markInternalIndel() { mInternalIndel = true; }

    public boolean discordantGroup() { return mDiscordantGroup; }
    public void markDiscordantGroup() { mDiscordantGroup = true; }

    public void setDepth(int depth) { mDepth = depth; }
    public int depth() { return mDepth; }

    public void addReadType(final ReadRecord read, final ReadType type)
    {
        if(ReadTypeReads.containsKey(type))
        {
            ReadTypeReads.get(type).add(read);
        }
    }

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
        return format("loc(%d:%d) frags(junc=%d exact=%d supp=%d) remotes(%d)",
                Position, Orientation, junctionFragmentCount(), exactSupportFragmentCount(), supportingFragmentCount(),
                RemoteJunctions.size());
    }
}
