package com.hartwig.hmftools.svtools.sv_prep;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SoftClipSide;

public class SvBucket
{
    private final int mId;

    private final List<ReadGroup> mReadGroups;
    private final List<ReadRecord> mSupportingReads; // for now only store as read

    private final List<JunctionData> mJunctions;
    private int mInitialSupportingReadCount;

    public SvBucket(final int id)
    {
        mId = id;
        mReadGroups = Lists.newArrayList();
        mSupportingReads = Lists.newArrayList();
        mJunctions = Lists.newArrayList();
        mInitialSupportingReadCount = 0;
    }

    public int id() { return mId; }

    public List<ReadGroup> readGroups() { return mReadGroups; }
    public List<ReadRecord> supportingReads() { return mSupportingReads; }
    public List<JunctionData> junctionPositions() { return mJunctions; }
    public int initialSupportingReadCount() { return mInitialSupportingReadCount; }

    public void addReadGroup(final ReadGroup readGroup)
    {
        mReadGroups.add(readGroup);
    }

    public void addSupportingRead(final ReadRecord read)
    {
        mSupportingReads.add(read);
    }

    public void setJunctionPositions()
    {
        for(ReadGroup readGroup : mReadGroups)
        {
            for(ReadRecord read : readGroup.reads())
            {
                SoftClipSide scSide = SoftClipSide.fromCigar(read.Cigar);

                if(scSide != null)
                {
                    if(scSide.isLeft())
                        addOrUpdateJunction(read.start(), NEG_ORIENT);
                    else
                        addOrUpdateJunction(read.end(), POS_ORIENT);
                }
            }
        }
    }

    private void addOrUpdateJunction(final int position, final byte orientation)
    {
        int index = 0;

        // use a binary search if count exeeds some level, but with buckets < 1000 bases likely unnecessary
        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionData.Position == position && junctionData.Orientation == orientation)
            {
                ++junctionData.ExactSupport;
                return;
            }
            else if(junctionData.Position >= position)
            {
                break;
            }

            ++index;
        }

        mJunctions.add(index, new JunctionData(position, orientation));
    }

    public void selectSupportingReads()
    {
        mInitialSupportingReadCount = mSupportingReads.size();

        List<ReadRecord> removeReads = Lists.newArrayList();
        int supplementaries = 0;

        for(ReadRecord read : mSupportingReads)
        {
            if(!readSupportsJunction(read))
            {
                removeReads.add(read);
            }
            else
            {
                // check for a 4th supp in an existing group
                ReadGroup readGroup = mReadGroups.stream().filter(x -> x.id().equals(read.Id)).findFirst().orElse(null);
                if(readGroup != null)
                {
                    readGroup.addRead(read);
                    removeReads.add(read);
                    ++supplementaries;
                }
            }
        }

        mInitialSupportingReadCount -= supplementaries;

        removeReads.forEach(x -> mSupportingReads.remove(x));
    }

    private boolean readSupportsJunction(final ReadRecord read)
    {
        int leftSoftClipPosition = read.Cigar.isLeftClipped() ? read.start() : 0;
        int rightSoftClipPosition = read.Cigar.isRightClipped() ? read.end() : 0;

        for(JunctionData junctionData : mJunctions)
        {
            // any length soft clipping at the same base
            if(junctionData.Position == leftSoftClipPosition && junctionData.Orientation == NEG_ORIENT)
            {
                ++junctionData.CandidateSupport;
                return true;
            }
            else if(junctionData.Position == rightSoftClipPosition && junctionData.Orientation == POS_ORIENT)
            {
                ++junctionData.CandidateSupport;
                return true;
            }

            // within 50 bases with exact sequence match in between the soft clip locations

        }

        return false;
    }
}
