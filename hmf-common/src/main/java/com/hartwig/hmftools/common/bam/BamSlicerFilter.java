package com.hartwig.hmftools.common.bam;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class BamSlicerFilter
{
    private final int mMinMappingQuality;
    private final boolean mKeepDuplicates;
    private final boolean mKeepSupplementaries;
    private final boolean mKeepSecondaries;
    private boolean mKeepHardClippedSecondaries; // when the -M was historically used in the BWA call
    private boolean mKeepUnmapped;

    public BamSlicerFilter(int minMappingQuality, boolean keepDuplicates, boolean keepSupplementaries, boolean keepSecondaries)
    {
        mMinMappingQuality = minMappingQuality;
        mKeepDuplicates = keepDuplicates;
        mKeepSupplementaries = keepSupplementaries;
        mKeepSecondaries = keepSecondaries;
        mKeepUnmapped = false;
        mKeepHardClippedSecondaries = false;
    }

    public void setKeepUnmapped()
    {
        mKeepUnmapped = true;
    }

    public void setKeepHardClippedSecondaries()
    {
        mKeepHardClippedSecondaries = true;
    }

    public boolean passesFilters(final SAMRecord record)
    {
        if(record.getMappingQuality() < mMinMappingQuality)
        {
            return false;
        }

        if(record.getReadUnmappedFlag() && !mKeepUnmapped)
        {
            return false;
        }

        if(record.isSecondaryAlignment() && !mKeepSecondaries)
        {
            if(!mKeepHardClippedSecondaries || !record.getCigar().containsOperator(CigarOperator.H))
            {
                return false;
            }
        }

        if(record.getSupplementaryAlignmentFlag() && !mKeepSupplementaries)
        {
            return false;
        }

        return !record.getDuplicateReadFlag() || mKeepDuplicates;
    }

    public boolean keepUnmapped()
    {
        return mKeepUnmapped;
    }
}
