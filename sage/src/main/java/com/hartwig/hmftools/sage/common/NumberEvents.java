package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public final class NumberEvents
{
    public static int calc(final SAMRecord record, final RefSequence refSequence)
    {
        int nm = rawNM(record, refSequence);

        int additionalIndels = 0;

        for(CigarElement cigarElement : record.getCigar())
        {
            switch(cigarElement.getOperator())
            {
                case D:
                case I:
                    additionalIndels += cigarElement.getLength() - 1;
                    break;
            }
        }

        return nm - additionalIndels;
    }

    public static double calcSoftClipAdjustment(final SAMRecord record)
    {
        return calcSoftClipAdjustment(leftSoftClipLength(record) + rightSoftClipLength(record));
    }

    public static double calcSoftClipAdjustment(int softClipLength)
    {
        return softClipLength > 0 ? max(1, softClipLength / SC_READ_EVENTS_FACTOR) : 0;
    }

    public static int rawNM(final SAMRecord record, final RefSequence refGenome)
    {
        Object nm = record.getAttribute(NUM_MUTATONS_ATTRIBUTE);
        if(nm instanceof Integer)
        {
            return (int) nm;
        }

        int offset = refGenome.Start - 1;
        return SequenceUtil.calculateSamNmTag(record, refGenome.Bases, offset);
    }

    public static int calcWithMnvRaw(int numberOfEvents, final String ref, final String alt)
    {
        // Number of events includes each SNV as an additional event. This unfairly penalises MNVs.
        int differentBases = 0;
        for(int i = 0; i < alt.length(); i++)
        {
            if(alt.charAt(i) != ref.charAt(i))
            {
                differentBases++;
            }
        }

        // We subtract one later when we actually use this value so we need to add one back in here to be consistent with SNVs and INDELs
        return numberOfEvents - differentBases + 1;
    }
}
