package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarBaseLength;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.alignedOrClipped;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public enum ReadValidReason
{
    OK,
    BASE_LENGTH,
    CIGAR_LENGTH,
    CIGAR_ELEMENTS,
    INVALID_BASE;

    public static ReadValidReason isValidRead(final SAMRecord record)
    {
        int baseLength = record.getReadBases().length;

        if(record.getBaseQualities().length != baseLength)
            return ReadValidReason.BASE_LENGTH;

        if(cigarBaseLength(record.getCigar()) != baseLength)
            return ReadValidReason.CIGAR_LENGTH;

        int cigarCount = record.getCigar().getCigarElements().size();
        if(cigarCount > 1)
        {
            for(int i = 0; i < record.getCigar().getCigarElements().size(); ++i)
            {
                CigarElement element = record.getCigar().getCigarElements().get(i);

                if(i == 0 || i == cigarCount - 1)
                {
                    if(!alignedOrClipped(element.getOperator()) && element.getOperator() != I)
                        return ReadValidReason.CIGAR_ELEMENTS;
                }
                else if(element.getOperator().isClipping())
                {
                    return ReadValidReason.CIGAR_ELEMENTS;
                }
            }
        }
        else if(record.getCigar().getCigarElements().get(0).getOperator() != M)
        {
            return ReadValidReason.CIGAR_ELEMENTS;
        }

        if(!hasValidBases(record))
            return ReadValidReason.INVALID_BASE;

        return ReadValidReason.OK;
    }

    public static boolean hasValidBases(final SAMRecord record)
    {
        for(int i = 0; i < record.getReadBases().length; ++i)
        {
            if(record.getReadBases()[i] == NO_BASE)
                return false;
        }

        return true;
    }
}
