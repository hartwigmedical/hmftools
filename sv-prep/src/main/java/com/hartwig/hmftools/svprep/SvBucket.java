package com.hartwig.hmftools.svprep;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvConstants.JUNCTION_SUPPORT_CAP;
import static com.hartwig.hmftools.svprep.SvConstants.LOW_BASE_QUALITY;
import static com.hartwig.hmftools.svprep.SvConstants.SUPPORTING_READ_DISTANCE;

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

    public boolean filterJunctions(int minJunctionSupport)
    {
        if(mReadGroups.isEmpty())
            return false;

        for(ReadGroup readGroup : mReadGroups)
        {
            for(ReadRecord read : readGroup.reads())
            {
                SoftClipSide scSide = SoftClipSide.fromCigar(read.mCigar);

                if(scSide != null)
                {
                    if(scSide.isLeft())
                        addOrUpdateJunction(read, NEG_ORIENT);
                    else
                        addOrUpdateJunction(read, POS_ORIENT);
                }
            }
        }

        selectSupportingReads();

        int index = 0;
        while(index < mJunctions.size())
        {
            if(mJunctions.get(index).totalSupport() < minJunctionSupport)
                mJunctions.remove(index);
            else
                ++index;
        }

        return !mJunctions.isEmpty();
    }

    private void addOrUpdateJunction(final ReadRecord read, final byte orientation)
    {
        int readPosition = orientation == NEG_ORIENT ? read.start() : read.end();

        int index = 0;

        // use a binary search if count exeeds some level, but with buckets < 1000 bases likely unnecessary
        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionData.Position == readPosition && junctionData.Orientation == orientation)
            {
                // should each read be tested for a match, and added?
                // junctionData.Reads.add(read);
                ++junctionData.ExactReads;
                return;
            }
            else if(junctionData.Position >= readPosition)
            {
                break;
            }

            ++index;
        }

        JunctionData junctionData = new JunctionData(readPosition, orientation, read);
        mJunctions.add(index, junctionData);
    }

    private void selectSupportingReads()
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
        boolean supportsJunction = false;

        if(mJunctions.size() < 20)
        {
            for(JunctionData junctionData : mJunctions)
            {
                // any length soft clipping at the same base
                if(supportsJunction(read, junctionData))
                {
                    ++junctionData.SupportReads;
                    supportsJunction = true;
                }
            }

            return supportsJunction;
        }

        int closeJunctionIndex = findJunctionIndex(read);

        if(closeJunctionIndex < 0)
            return false;

        // check up and down from this location
        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);

            int index = searchUp ? closeJunctionIndex + 1 : closeJunctionIndex - 1;

            while(index >= 0 && index < mJunctions.size())
            {
                JunctionData junctionData = mJunctions.get(index);

                if(!readWithinJunctionRange(read, junctionData))
                    break;

                if(junctionData.SupportReads < JUNCTION_SUPPORT_CAP) // TEMP to limit processing
                {
                    if(supportsJunction(read, junctionData))
                    {
                        ++junctionData.SupportReads;
                        supportsJunction = true;
                    }
                }

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }

        return supportsJunction;
    }

    private int findJunctionIndex(final ReadRecord read)
    {
        // binary search on junctions if the collection gets too large
        int currentIndex = mJunctions.size() / 2;
        int lowerIndex = 0;
        int upperIndex = mJunctions.size() - 1;

        while(true)
        {
            JunctionData junctionData = mJunctions.get(currentIndex);

            if(readWithinJunctionRange(read, junctionData))
                return currentIndex;

            if(positionWithin(read.start(), read.end(), junctionData.Position))
                return -1;

            if(read.end() < junctionData.Position)
            {
                // search lower
                if(lowerIndex + 1 == currentIndex)
                    return currentIndex;

                upperIndex = currentIndex;
                currentIndex = (lowerIndex + upperIndex) / 2;
            }
            else if(read.start() > junctionData.Position)
            {
                // search higher
                if(currentIndex + 1 == upperIndex)
                    return currentIndex;

                lowerIndex = currentIndex;
                currentIndex = (lowerIndex + upperIndex) / 2;
            }
            else
            {
                break;
            }
        }

        return -1;
    }

    private static boolean readWithinJunctionRange(final ReadRecord read, final JunctionData junctionData)
    {
        if(read.cigar().isRightClipped())
        {
            if(abs(read.end() - junctionData.Position) <= SUPPORTING_READ_DISTANCE)
                return true;
        }

        if(read.cigar().isLeftClipped())
        {
            if(abs(read.start() - junctionData.Position) <= SUPPORTING_READ_DISTANCE)
                return true;
        }

        return false;
    }

    public static boolean supportsJunction(final ReadRecord read, final JunctionData junctionData)
    {
        /*
        eg if there is a candidate read with 101M50S at base 1000, how does a supporting read need to align and match?
         In this case the soft clip is at 1050. Check all reads which have a soft clip on the same side within +/- 50 bases.
        eg. you may find a read with soft clip at 1047 with 10S.
        Check if the 10 bases match the last 3 M and first 7S of the original read allowing for low base qual mismatches.
         If they do then that counts as an additional read even though it is much shorter soft clip and does not exactly match the base.
         */
        final ReadRecord juncRead = junctionData.InitialRead;

        if(junctionData.Orientation == POS_ORIENT)
        {
            if(!read.cigar().isRightClipped())
                return false;

            int readRightPos = read.end();

            if(readRightPos == junctionData.Position)
                return true;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readRightPos - junctionData.Position) > SUPPORTING_READ_DISTANCE)
                return false;

            int scLength = read.cigar().getLastCigarElement().getLength();
            int readLength = read.readBases().length();
            int readEndPosIndex = readLength - scLength - 1;

            int juncReadLength = juncRead.readBases().length();
            int juncReadScLength = juncRead.cigar().getLastCigarElement().getLength();
            int juncReadEndPosIndex = juncReadLength - juncReadScLength - 1;
            int endPosDiff = juncRead.end() - readRightPos;

            int junctionReadOffset = juncReadEndPosIndex - readEndPosIndex - endPosDiff;

            // test against all the read's right soft-clipped bases
            for(int i = readLength - scLength; i < readLength; ++i)
            {
                char readBase = read.readBases().charAt(i);

                int juncIndex = i + junctionReadOffset;
                if(juncIndex < 0 || juncIndex >= juncReadLength)
                    return false;

                char juncReadBase = juncRead.readBases().charAt(juncIndex);

                if(readBase == juncReadBase)
                    continue;

                if(read.baseQualities()[i] < LOW_BASE_QUALITY || juncRead.baseQualities()[juncIndex] < LOW_BASE_QUALITY)
                    continue;

                return false;
            }
        }
        else
        {
            if(!read.cigar().isLeftClipped())
                return false;

            int readLeftPos = read.start();

            if(readLeftPos == junctionData.Position)
                return true;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readLeftPos - junctionData.Position) > SUPPORTING_READ_DISTANCE)
                return false;

            // test for a base match for the read's soft-clipped bases, allow for low-qual matches

            // read: SC length -> start position
            // junc: SC length -> start position
            // junc read index = sc length diff - position diff

            int scLength = read.cigar().getFirstCigarElement().getLength();

            int juncReadScLength = juncRead.cigar().getFirstCigarElement().getLength();
            int posOffset = juncRead.start() - readLeftPos;
            int softClipDiff = juncReadScLength - scLength;
            int junctionReadOffset = softClipDiff - posOffset;
            int juncReadLength = juncRead.readBases().length();

            for(int i = 0; i < scLength; ++i)
            {
                char readBase = read.readBases().charAt(i);

                int juncIndex = i + junctionReadOffset;
                if(juncIndex < 0 || juncIndex >= juncReadLength)
                    return false;

                char juncReadBase = juncRead.readBases().charAt(juncIndex);

                if(readBase == juncReadBase)
                    continue;

                if(read.baseQualities()[i] < LOW_BASE_QUALITY || juncRead.baseQualities()[juncIndex] < LOW_BASE_QUALITY)
                    continue;

                return false;
            }
        }

        return true;
    }
}
