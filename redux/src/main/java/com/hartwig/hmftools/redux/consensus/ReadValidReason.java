package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarBaseLength;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.alignedOrClipped;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.redux.ReduxConfig;

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

        if(!hasValidCigar(record.getCigar().getCigarElements()))
            return ReadValidReason.CIGAR_ELEMENTS;

        if(!hasValidBases(record))
            return ReadValidReason.INVALID_BASE;

        return ReadValidReason.OK;
    }

    public static boolean hasValidCigar(final List<CigarElement> cigarElements)
    {
        int cigarCount = cigarElements.size();

        if(cigarCount > 1)
        {
            for(int i = 0; i < cigarCount; ++i)
            {
                CigarElement element = cigarElements.get(i);

                if(i == 0 || i == cigarCount - 1)
                {
                    if(!alignedOrClipped(element.getOperator()) && element.getOperator() != I)
                        return false;
                }
                else if(element.getOperator().isClipping())
                {
                    return false;
                }
            }
        }
        else if(cigarElements.get(0).getOperator() != M)
        {
            return false;
        }

        return true;
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

    public static List<Integer> getInvalidBases(final SAMRecord record)
    {
        List<Integer> invalidBases = Lists.newArrayList();

        for(int i = 0; i < record.getReadBases().length; ++i)
        {
            if(record.getReadBases()[i] == NO_BASE)
                invalidBases.add(i);
        }

        return invalidBases;
    }

    public static void checkIsValidRead(final SAMRecord read)
    {
        ReadValidReason validReason = isValidRead(read);

        if(validReason != ReadValidReason.OK)
        {
            RD_LOGGER.error("invalid read({}): {}", validReason, readToString(read));
        }
    }
}
