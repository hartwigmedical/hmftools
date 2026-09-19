package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class MappedCoords
{
    private final int mOriginalCount;
    private final List<BaseRegion> mAlignments;

    private boolean[] mInferredAlignmentAdded;
    private int[] mSoftClipRegionsMatched;

    public MappedCoords(final List<BaseRegion> alignments)
    {
        mOriginalCount = alignments.size();
        mAlignments = alignments;
        mInferredAlignmentAdded = null;
        mSoftClipRegionsMatched = null;
    }

    public int alignmentCount() { return mAlignments.size(); }
    public int originalAlignmentCount() { return mOriginalCount; }

    public List<BaseRegion> alignments() { return mAlignments; }

    public List<BaseRegion> alignmentsWithoutInferred()
    {
        if(mInferredAlignmentAdded == null)
            return mAlignments;

        List<BaseRegion> regions = Lists.newArrayListWithCapacity(mOriginalCount);

        int startIndex = mInferredAlignmentAdded[SE_START] ? 1 : 0;

        int endIndex = startIndex + mOriginalCount - 1;
        for(int i = startIndex; i <= endIndex; ++i)
        {
            regions.add(mAlignments.get(i));
        }

        return regions;
    }

    public BaseRegion lowestAlignment(boolean includeInferred)
    {
        if(includeInferred || mInferredAlignmentAdded == null || !mInferredAlignmentAdded[SE_START])
            return mAlignments.get(0);

        return mAlignments.get(1);
    }

    public BaseRegion highestAlignment(boolean includeInferred)
    {
        if(includeInferred || mInferredAlignmentAdded == null || !mInferredAlignmentAdded[SE_END])
            return mAlignments.get(mAlignments.size() - 1);

        return mAlignments.get(mAlignments.size() - 2);
    }

    public BaseRegion regionByIndex(int index) { return mAlignments.get(index); }

    public static final int INVALID_INDEX = -1;

    public int findRegionIndex(final RegionReadData region)
    {
        for(int i = 0; i < mAlignments.size(); ++i)
        {
            BaseRegion readSection = mAlignments.get(i);

            if(positionsOverlap(readSection.start(), readSection.end(), region.start(), region.end()))
                return i;
        }

        return INVALID_INDEX;
    }

    public boolean lowerInferredAlignmentAdded() { return mInferredAlignmentAdded != null ? mInferredAlignmentAdded[SE_START] : false; }
    public boolean upperInferredAlignmentAdded() { return mInferredAlignmentAdded != null ? mInferredAlignmentAdded[SE_END] : false; }
    public boolean inferredAlignmentAdded(int seIndex) { return mInferredAlignmentAdded != null ? mInferredAlignmentAdded[seIndex] : false; }

    public void addInferredRegion(boolean isLower, int posStart, int posEnd)
    {
        if(mInferredAlignmentAdded == null)
            mInferredAlignmentAdded = new boolean[] {false, false};

        if(isLower)
        {
            if(!mInferredAlignmentAdded[SE_START])
            {
                mInferredAlignmentAdded[SE_START] = true;
                mAlignments.add(0, new BaseRegion(posStart, posEnd));
            }
            else
            {
                // lengthen the new region if required
                BaseRegion existing = mAlignments.get(0);
                existing.setStart(min(existing.start(), posStart));
            }
        }
        else
        {
            if(!mInferredAlignmentAdded[SE_END])
            {
                mInferredAlignmentAdded[SE_END] = true;
                mAlignments.add(new BaseRegion(posStart, posEnd));
            }
            else
            {
                BaseRegion existing = mAlignments.get(mAlignments.size() - 1);
                existing.setEnd(max(existing.end(), posEnd));
            }
        }
    }

    public boolean alignmentsOverlap(int posStart, int posEnd)
    {
        return mAlignments.stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x.start(), x.end()));
    }

    public boolean alignmentsWithin(int posStart, int posEnd)
    {
        return mAlignments.stream().anyMatch(x -> positionsWithin(x.start(), x.end(), posStart, posEnd));
    }

    public int getCoordsBoundary(int se)
    {
        return se == SE_START ? mAlignments.get(0).start() : mAlignments.get(mAlignments.size() - 1).end();
    }

    public void addSoftClipRegionMatched(boolean onStart, int count)
    {
        if(mSoftClipRegionsMatched == null)
        {
            mSoftClipRegionsMatched = new int[] {0, 0};
        }

        mSoftClipRegionsMatched[onStart ? SE_START : SE_END] += count;
    }

    public boolean isSoftClipRegionMatched(int seIndex)
    {
        return mSoftClipRegionsMatched != null ? mSoftClipRegionsMatched[seIndex] > 0 : false;
    }

    public int softClipRegionsMatched(int seIndex)
    {
        return mSoftClipRegionsMatched != null ? mSoftClipRegionsMatched[seIndex] : 0;
    }

    public String toString()
    {
        if(mAlignments.size() == 1)
            return mAlignments.get(0).toString();

        String alignmentsStr = mAlignments.stream().map(x -> x.toString()).collect(Collectors.joining(";"));

        if(mSoftClipRegionsMatched == null && mInferredAlignmentAdded == null)
            return format("%s", alignmentsStr);

        return format("%s inferred(lower=%s upper=%d) regionsMatched(lower=%d upper=%d)",
                alignmentsStr, mInferredAlignmentAdded[SE_START], mInferredAlignmentAdded[SE_END],
                mSoftClipRegionsMatched[SE_START], mSoftClipRegionsMatched[SE_END]);
    }

    public static MappedCoords build(final List<CigarElement> cigarElements, int posStart)
    {
        int splitCount = (int)cigarElements.stream().filter(x -> x.getOperator() == CigarOperator.N).count();

        List<BaseRegion> alignments = Lists.newArrayListWithCapacity(1 * splitCount);

        BaseRegion currentRegion = null;

        int posOffset = 0;
        boolean continueRegion = false;

        for(CigarElement element : cigarElements)
        {
            if(element.getOperator() == S)
                continue;

            if(element.getOperator() == I)
            {
                continueRegion = true;
                continue;
            }

            if(element.getOperator() == D)
            {
                // don't break alignments for deletes
                posOffset += element.getLength();
                continueRegion = true;
            }
            else if(element.getOperator() == CigarOperator.N)
            {
                // break on splits
                posOffset += element.getLength();
                continueRegion = false;
            }
            else if(element.getOperator() == CigarOperator.M)
            {
                int readStartPos = posStart + posOffset;
                int readEndPos = readStartPos + element.getLength() - 1;

                if(continueRegion && currentRegion != null)
                {
                    currentRegion.setEnd(readEndPos);
                }
                else
                {
                    currentRegion = new BaseRegion(readStartPos, readEndPos);
                    alignments.add(currentRegion);
                }

                posOffset += element.getLength();
                continueRegion = false;
            }
        }

        return new MappedCoords(alignments);
    }
}
