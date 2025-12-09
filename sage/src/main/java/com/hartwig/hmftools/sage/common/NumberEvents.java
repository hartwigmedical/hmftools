package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;

import static htsjdk.samtools.util.FileExtensions.CRAM;

import java.util.List;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public final class NumberEvents
{
    public static boolean RECOMPUTE_MISSING_NM = false;

    public static void setRecomputeNumMutations(final List<String> bamFiles)
    {
        // assumption is that CRAMs have dropped NM for compression and if missing it needs to be recomputed
        if(bamFiles.stream().anyMatch(x -> x.endsWith(CRAM)))
            RECOMPUTE_MISSING_NM = true;
    }

    public static int calcAdjustedNumMutations(final SAMRecord record, final RefSequence refSequence)
    {
        int nm = getOrCalcNm(record, refSequence);

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

    private static int getOrCalcNm(final SAMRecord record, final RefSequence refSequence)
    {
        Object nm = record.getAttribute(NUM_MUTATONS_ATTRIBUTE);

        if(nm == null && !RECOMPUTE_MISSING_NM)
            return 0;

        if(nm instanceof Integer)
            return (int) nm;

        int offset = refSequence.Start - 1;
        return SequenceUtil.calculateSamNmTag(record, refSequence.Bases, offset);
    }

    public static double calcSoftClipAdjustment(final SAMRecord record)
    {
        return calcSoftClipAdjustment(leftSoftClipLength(record) + rightSoftClipLength(record));
    }

    public static double calcSoftClipAdjustment(int softClipLength)
    {
        return softClipLength > 0 ? max(1, softClipLength / SC_READ_EVENTS_FACTOR) : 0;
    }

    public static int calcWithMnvRaw(int numberOfEvents, final String ref, final String alt)
    {
        // number of events includes each SNV as an additional event - this unfairly penalises MNVs
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
