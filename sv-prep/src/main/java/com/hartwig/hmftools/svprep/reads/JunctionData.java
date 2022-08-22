package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.getSoftClippedBases;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class JunctionData
{
    public final int Position;
    public final byte Orientation;

    public final List<ReadGroup> JunctionGroups; // with a read matching the junction
    public final List<ReadGroup> SupportingGroups;
    public final List<ReadGroup> ExactSupportGroups;
    public final List<RemoteJunction> RemoteJunctions;

    public final Map<ReadType,List<ReadRecord>> ReadTypeReads;

    private ReadRecord mTopJunctionRead;
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

        mTopJunctionRead = read; // initially set to the first

        mHotspot = false;
        mInternalIndel = false;
        mDiscordantGroup = false;
        mDepth = 0;
    }

    public boolean isExisting() { return mTopJunctionRead == null; }

    public ReadRecord topJunctionRead() { return mTopJunctionRead; }

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

    public void setInitialRead(int minSoftClipHighQual)
    {
        if(mInternalIndel)
            return;

        // select the junction read with the highest number of high-qual bases in the soft-clip
        int maxHighQualBases = 0;
        ReadRecord topRead = null;

        boolean useLeftSoftClip = Orientation == NEG_ORIENT;

        List<ReadRecord> junctionReads = ReadTypeReads.get(ReadType.JUNCTION);

        if(junctionReads.size() == 1)
        {
            mTopJunctionRead = junctionReads.get(0);
            return;
        }

        for(ReadRecord read : junctionReads)
        {
            String scBases = getSoftClippedBases(read.record(), useLeftSoftClip);

            if(scBases.length() < maxHighQualBases)
                continue;

            final byte[] baseQualities = read.record().getBaseQualities();
            int scRangeStart = useLeftSoftClip ? 0 : baseQualities.length - scBases.length();
            int scRangeEnd = useLeftSoftClip ? scBases.length() : baseQualities.length;

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
        return format("loc(%d:%d) frags(junc=%d exact=%d supp=%d) remotes(%d)",
                Position, Orientation, junctionFragmentCount(), exactSupportFragmentCount(), supportingFragmentCount(),
                RemoteJunctions.size());
    }
}
