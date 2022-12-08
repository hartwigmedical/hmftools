package com.hartwig.hmftools.bamtools.markdups;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.ReadGroup;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroup
{
    private final ChrBaseRegion mLowerReadRange;
    private final ChrBaseRegion mUpperReadRange;
    private final List<ReadGroup> mReadGroups;

    public DuplicateGroup(final ChrBaseRegion lowerCoords, final ChrBaseRegion upperCoords, final ReadGroup first, final ReadGroup second)
    {
        mReadGroups = Lists.newArrayList();
        mReadGroups.add(first);
        mReadGroups.add(second);

        mLowerReadRange = lowerCoords;
        mUpperReadRange = upperCoords;
    }

    public ChrBaseRegion lowerCoords() { return mLowerReadRange; }
    public ChrBaseRegion upperCoords() { return mUpperReadRange; }
    public List<ReadGroup> readGroups() { return mReadGroups; }

    public static ChrBaseRegion lowerCoords(final ReadGroup readGroup)
    {
        SAMRecord lowerRead = readGroup.reads().get(0);
        return new ChrBaseRegion(lowerRead.getContig(), lowerRead.getAlignmentStart(), lowerRead.getAlignmentEnd());
    }

    public static ChrBaseRegion upperCoords(final ReadGroup readGroup)
    {
        SAMRecord lowerRead = readGroup.reads().get(0);
        if(readGroup.reads().size() > 1)
        {
            SAMRecord upperRead = readGroup.reads().get(1);
            return new ChrBaseRegion(upperRead.getContig(), upperRead.getAlignmentStart(), upperRead.getAlignmentEnd());
        }
        else if(lowerRead.getMateUnmappedFlag())
        {
            return null;
        }
        else
        {
            return new ChrBaseRegion(lowerRead.getMateReferenceName(), lowerRead.getMateAlignmentStart(), 0);
        }
    }

    public static boolean upperCoordsMatch(final ChrBaseRegion first, final ChrBaseRegion second)
    {
        if(first == null && second == null)
            return true;

        if(first == null || second == null)
            return false;

        return first.matches(second);
    }
}
