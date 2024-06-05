package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class JunctionData
{
    public final int Position;
    public final Orientation Orient;

    private final boolean mIsExisting;
    public final List<ReadGroup> JunctionGroups; // with a read matching the junction
    public final List<ReadGroup> SupportingGroups;
    public final List<ReadGroup> ExactSupportGroups;
    public final List<RemoteJunction> RemoteJunctions;

    public final Map<ReadType,List<PrepRead>> ReadTypeReads;

    private PrepRead mTopJunctionRead;
    private boolean mInternalIndel;
    private boolean mDiscordantGroup;
    private boolean mHotspot;

    public JunctionData(final int position, final Orientation orientation, final PrepRead read)
    {
        Position = position;
        Orient = orientation;

        mIsExisting = read == null;
        JunctionGroups = Lists.newArrayList();
        SupportingGroups = Lists.newArrayList();
        ExactSupportGroups = Lists.newArrayList();
        RemoteJunctions = Lists.newArrayList();
        ReadTypeReads = Maps.newHashMap();

        ReadTypeReads.put(ReadType.JUNCTION, Lists.newArrayList());
        ReadTypeReads.put(ReadType.SUPPORT, Lists.newArrayList());
        ReadTypeReads.put(ReadType.EXACT_SUPPORT, Lists.newArrayList());

        mTopJunctionRead = read; // initially set to the first

        mHotspot = false;
        mInternalIndel = false;
        mDiscordantGroup = false;
    }

    public boolean isExisting() { return mIsExisting; }
    public boolean isForward() { return Orient.isForward(); }
    public boolean isReverse() { return Orient.isReverse(); }

    public PrepRead topJunctionRead() { return mTopJunctionRead; }

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

    public void addReadType(final PrepRead read, final ReadType type)
    {
        if(ReadTypeReads.containsKey(type))
        {
            ReadTypeReads.get(type).add(read);
        }
    }

    public void setInitialRead(int minSoftClipHighQual)
    {
        if(mInternalIndel)
            return;

        // select the junction read with the highest number of high-qual bases in the soft-clip
        int maxHighQualBases = 0;
        PrepRead topRead = null;

        boolean useLeftSoftClip = Orient.isReverse();

        List<PrepRead> junctionReads = ReadTypeReads.get(ReadType.JUNCTION);

        if(junctionReads.size() == 1)
        {
            mTopJunctionRead = junctionReads.get(0);
            return;
        }

        for(PrepRead read : junctionReads)
        {
            int scLength = useLeftSoftClip ? read.leftClipLength() : read.rightClipLength();

            if(scLength < maxHighQualBases)
                continue;

            final byte[] baseQualities = read.record().getBaseQualities();
            int scRangeStart = useLeftSoftClip ? 0 : baseQualities.length - scLength;
            int scRangeEnd = useLeftSoftClip ? scLength : baseQualities.length;

            int aboveQual = 0;
            for(int i = scRangeStart; i < scRangeEnd; ++i)
            {
                if(baseQualities[i] >= minSoftClipHighQual)
                    ++aboveQual;
            }

            if(aboveQual > maxHighQualBases)
            {
                topRead = read;
                maxHighQualBases = aboveQual;
            }
        }

        if(topRead != null)
            mTopJunctionRead = topRead;
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
        return format("loc(%d:%d) frags(junc=%d exact=%d supp=%d) remotes(%d) disc(%s) indel(%s)",
                Position, Orient.asByte(), junctionFragmentCount(), exactSupportFragmentCount(), supportingFragmentCount(),
                RemoteJunctions.size(), mDiscordantGroup, mInternalIndel);
    }
}
