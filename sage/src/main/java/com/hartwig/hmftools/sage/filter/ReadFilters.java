package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.sage.SageConstants.CHIMERIC_FRAGMENT_LENGTH_MAX;

import htsjdk.samtools.SAMRecord;

public final class ReadFilters
{
    public static boolean isChimericRead(final SAMRecord record)
    {
        if(record.getReadPairedFlag())
        {
            if(record.getMateUnmappedFlag())
                return true;

            // inter-chromosomal
            if(!record.getReferenceName().equals(record.getMateReferenceName()))
                return true;

            // inversion
            if(record.getReadNegativeStrandFlag() == mateNegativeStrand(record))
                return true;
        }

        // or a fragment length outside the expected maximum
        if(abs(record.getInferredInsertSize()) > CHIMERIC_FRAGMENT_LENGTH_MAX)
            return true;

        return false;
    }
}
