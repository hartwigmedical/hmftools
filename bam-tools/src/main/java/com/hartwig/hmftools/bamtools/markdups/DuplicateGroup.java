package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.round;

import java.util.List;
import java.util.stream.Collectors;

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

    public ReadGroup findPrimaryGroup()
    {
        List<ReadGroup> nonDupGroups = mReadGroups.stream().filter(x -> !hasDuplicates(x)).collect(Collectors.toList());

        if(nonDupGroups.size() == 1)
            return nonDupGroups.get(0);

        ReadGroup maxGroup = null;
        int maxBaseQual = 0;

        for(ReadGroup readGroup : mReadGroups)
        {
            int groupBaseQual = calcBaseQualTotal(readGroup);

            if(groupBaseQual > maxBaseQual)
            {
                maxBaseQual = groupBaseQual;
                maxGroup = readGroup;
            }
        }

        return maxGroup;
    }

    public static int calcBaseQualTotal(final ReadGroup readGroup)
    {
        int readBaseCount = 0;
        int readBaseQualTotal = 0;

        for(SAMRecord read : readGroup.reads())
        {
            if(read.getSupplementaryAlignmentFlag())
                continue;

            for(int i = 0; i < read.getBaseQualities().length; ++i)
            {
                ++readBaseCount;
                readBaseQualTotal += read.getBaseQualities()[i];
            }
        }

        return readBaseCount > 0 ? (int)round(readBaseQualTotal / (double)readBaseCount) : 0;
    }

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

    public static boolean hasDuplicates(final ReadGroup readGroup)
    {
        return readGroup.reads().stream().anyMatch(x -> x.getDuplicateReadFlag());
    }
}
