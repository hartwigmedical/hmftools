package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.aboveLowBaseQual;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class JunctionData
{
    public final int Position;
    public final Orientation Orient;

    private List<ReadGroup> mJunctionGroups; // with a read matching the junction
    private List<ReadGroup> mSupportingGroups;
    private List<ReadGroup> mExactSupportGroups;
    private List<RemoteJunction> mRemoteJunctions;

    private final Map<ReadType,List<PrepRead>> mReadTypeReads;

    private PrepRead mTopJunctionRead;
    private boolean mInternalIndel;
    private boolean mDiscordantGroup;
    private boolean mHotspot;
    private int mDepth;

    private JunctionData mLinkedIndel;

    public JunctionData(final int position, final Orientation orientation, final PrepRead read)
    {
        Position = position;
        Orient = orientation;

        mJunctionGroups = null;
        mSupportingGroups = null;
        mExactSupportGroups = null;
        mRemoteJunctions = null;
        mReadTypeReads = Maps.newHashMap();

        mReadTypeReads.put(ReadType.JUNCTION, Lists.newArrayList());
        mReadTypeReads.put(ReadType.SUPPORT, Lists.newArrayList());
        mReadTypeReads.put(ReadType.EXACT_SUPPORT, Lists.newArrayList());

        mTopJunctionRead = read; // initially set to the first

        mHotspot = false;
        mInternalIndel = false;
        mDiscordantGroup = false;
        mDepth = 0;
        mLinkedIndel = null;
    }

    public boolean isForward() { return Orient.isForward(); }
    public boolean isReverse() { return Orient.isReverse(); }

    public PrepRead topJunctionRead() { return mTopJunctionRead; }

    public void addJunctionReadGroup(final ReadGroup readGroup)
    {
        if(mJunctionGroups == null)
            mJunctionGroups = Lists.newArrayList();

        mJunctionGroups.add(readGroup);
    }

    public void addExactSupportGroup(final ReadGroup readGroup)
    {
        if(mExactSupportGroups == null)
            mExactSupportGroups = Lists.newArrayList();

        mExactSupportGroups.add(readGroup);
    }

    public void addSupportingGroup(final ReadGroup readGroup)
    {
        if(mSupportingGroups == null)
            mSupportingGroups = Lists.newArrayList();

        mSupportingGroups.add(readGroup);
    }

    public List<ReadGroup> junctionGroups() { return mJunctionGroups != null ? mJunctionGroups : Collections.emptyList(); }
    public List<ReadGroup> supportingGroups() { return mSupportingGroups != null ? mSupportingGroups : Collections.emptyList(); }
    public List<ReadGroup> exactSupportGroups() { return mExactSupportGroups != null ? mExactSupportGroups : Collections.emptyList(); }

    public Map<ReadType,List<PrepRead>> readTypeReads() { return mReadTypeReads; }

    public int junctionFragmentCount() { return mJunctionGroups != null ? mJunctionGroups.size() : 0; }
    public int supportingFragmentCount() { return mSupportingGroups != null ? mSupportingGroups.size() : 0; }
    public int exactSupportFragmentCount() { return mExactSupportGroups != null ? mExactSupportGroups.size() : 0; }

    public List<RemoteJunction> remoteJunctions() { return mRemoteJunctions != null ? mRemoteJunctions : Collections.emptyList(); }

    public boolean hotspot() { return mHotspot; }
    public void markHotspot() { mHotspot = true; }

    public boolean internalIndel() { return mInternalIndel; }
    public void markInternalIndel() { mInternalIndel = true; }

    public boolean discordantGroup() { return mDiscordantGroup; }
    public void markDiscordantGroup() { mDiscordantGroup = true; }

    public int depth() { return mDepth; }
    public void setDepth(double depth) { mDepth = (int)depth; }

    public void addReadType(final PrepRead read, final ReadType type)
    {
        if(mReadTypeReads.containsKey(type))
        {
            mReadTypeReads.get(type).add(read);
        }
    }

    public void setInitialRead()
    {
        if(mInternalIndel)
            return;

        // select the junction read with the highest number of high-qual bases in the soft-clip
        int maxHighQualBases = 0;
        PrepRead topRead = null;

        boolean useLeftSoftClip = Orient.isReverse();

        List<PrepRead> junctionReads = mReadTypeReads.get(ReadType.JUNCTION);

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

            if(scRangeEnd - scRangeStart + 1 < maxHighQualBases)
                continue;

            int aboveQual = 0;
            for(int i = scRangeStart; i < scRangeEnd; ++i)
            {
                if(aboveLowBaseQual(baseQualities[i]))
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
        if(mRemoteJunctions == null)
            mRemoteJunctions = Lists.newArrayList();

        RemoteJunction matched = mRemoteJunctions.stream().filter(x -> x.matches(remoteJunction)).findFirst().orElse(null);
        if(matched != null)
        {
            ++matched.Fragments;
            return;
        }

        mRemoteJunctions.add(remoteJunction);
    }

    public JunctionData linkedIndel() { return mLinkedIndel; }
    public void setLinkedIndel(final JunctionData junctionData) { mLinkedIndel = junctionData; }
    public boolean hasLinkedIndel() { return mLinkedIndel != null; }

    public String toString()
    {
        return format("loc(%d:%d) frags(junc=%d exact=%d supp=%d) remotes(%d) disc(%s) indel(%s)",
                Position, Orient.asByte(), junctionFragmentCount(), exactSupportFragmentCount(), supportingFragmentCount(),
                mRemoteJunctions != null ? mRemoteJunctions.size() : 0, mDiscordantGroup, mInternalIndel);
    }
}
