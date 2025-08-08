package com.hartwig.hmftools.common.sequencing;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.PHRED_OFFSET;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;

import com.hartwig.hmftools.common.bam.ConsensusType;

import htsjdk.samtools.SAMRecord;

public final class UltimaBamUtils
{
    public static final byte ULTIMA_MAX_QUAL = 40;
    public static final byte ULTIMA_INVALID_QUAL = -1;

    public static final String ULTIMA_TP_TAG = "tp";
    public static final String ULTIMA_T0_TAG = "t0";

    private static final String PPM_STRAND_AS = "as";
    private static final String PPM_STRAND_TS = "ts";
    private static final String PPM_STRAND_AE = "ae";
    private static final String PPM_STRAND_TE = "te";

    private static final int PPM_STRAND_MIN_SUM = 4;
    private static final int PPM_STRAND_MAX_SUM = 8;
    private static final double PPM_STRAND_BALANCED_LOW = 0.27;
    private static final double PPM_STRAND_BALANCED_HIGH = 1 - PPM_STRAND_BALANCED_LOW;

    public static final byte[] CYCLE_BASES = new byte[] { (byte)'T', (byte)'G', (byte)'C', (byte)'A' };

    public static boolean isHighBaseQual(final byte qual)
    {
        return qual >= ULTIMA_MAX_QUAL;
    }

    public static byte[] extractTpValues(final SAMRecord record)
    {
        return record.getByteArrayAttribute(ULTIMA_TP_TAG);
    }

    public static byte[] extractT0Values(final SAMRecord record)
    {
        byte[] t0Values = record.getStringAttribute(ULTIMA_T0_TAG).getBytes();
        for(int i = 0; i < t0Values.length; ++i)
        {
            t0Values[i] -= PHRED_OFFSET;
        }
        return t0Values;
    }

    public static ConsensusType deriveConsensusType(final SAMRecord record)
    {
        //  example: as:i:3  ts:i:4  ae:i:3  te:i:3
        int as = getIntegerAttribute(record, PPM_STRAND_AS, 0);
        int ts = getIntegerAttribute(record, PPM_STRAND_TS, 0);
        int ae = getIntegerAttribute(record, PPM_STRAND_AE, 0);
        int te = getIntegerAttribute(record, PPM_STRAND_TE, 0);

        double sumStart = as + ts;
        double sumEnd = ae + te;

        if(sumStart < PPM_STRAND_MIN_SUM || sumStart > PPM_STRAND_MAX_SUM || sumEnd < PPM_STRAND_MIN_SUM || sumEnd > PPM_STRAND_MAX_SUM)
            return ConsensusType.NONE;

        double ratioStart = as / sumStart;
        double ratioEnd = ae / sumEnd;

        if(ratioStart < PPM_STRAND_BALANCED_LOW || ratioStart > PPM_STRAND_BALANCED_HIGH)
            return ConsensusType.NONE;

        if(ratioEnd < PPM_STRAND_BALANCED_LOW || ratioEnd > PPM_STRAND_BALANCED_HIGH)
            return ConsensusType.NONE;

        return ConsensusType.DUAL;
    }

    public static byte calcTpBaseQual(final SAMRecord record, int indexStart, int indexEnd, int tpSearchValue)
    {
        byte[] tpValues = extractTpValues(record);

        int qualValue1 = -1;
        int qualValue2 = -1;

        for(int i = indexStart; i <= min(indexEnd, tpValues.length - 1); ++i)
        {
            if(tpValues[i] == tpSearchValue)
            {
                if(qualValue1 < 0)
                {
                    qualValue1 = record.getBaseQualities()[i];
                }
                else
                {
                    qualValue2 = record.getBaseQualities()[i];
                    break;
                }
            }
        }

        if(qualValue1 < 0)
        {
            // if bases aren't found, use middle base(s) as: min(40, Q + 6 * abs(required_tp – T))
            int middleIndex = indexStart + (indexEnd - indexStart) / 2;
            qualValue1 = record.getBaseQualities()[middleIndex];
            byte tpValue = tpValues[middleIndex];

            return (byte)min(ULTIMA_MAX_QUAL, qualValue1 + 6 * abs(tpSearchValue - tpValue));
        }

        if(qualValue2 < 0)
            return (byte)qualValue1;

        // equivalent to adding their logs, ie P(combined) = 10^(-qual/10) + 10^(-qual/10), combined qual = -10 * log10(P(combined))
        return (byte)(qualValue1 - 3);
    }

    public static int cycleCount(final byte firstBase, final byte lastBase, final byte innerBase)
    {
        if(baseIndex(firstBase) == -1 || baseIndex(lastBase) == -1 || baseIndex(innerBase) == -1)
            return -1;  // any base is not TGCA
        byte[] bases = new byte[] {firstBase, innerBase, lastBase};
        int baseIndex = 0;
        int cycleIndex = 0;
        int cycleCount = 1;

        while(baseIndex < bases.length)
        {
            if(CYCLE_BASES[cycleIndex] == bases[baseIndex])
                baseIndex += 1;
            else
                cycleIndex += 1;

            if(cycleIndex == CYCLE_BASES.length)
            {
                cycleIndex = 0;
                cycleCount += 1;
            }
        }
        return cycleCount;
    }

    public static boolean isBaseInCycle(final byte startBase, final byte endBase, final byte testBase)
    {
        if(startBase == endBase || startBase == testBase)
            return false;

        boolean started = false;
        int index = 0;

        while(true)
        {
            if(!started)
            {
                if(CYCLE_BASES[index] == startBase)
                    started = true;
            }
            else
            {
                if(CYCLE_BASES[index] == endBase || CYCLE_BASES[index] == startBase)
                    return false;
                else if(CYCLE_BASES[index] == testBase)
                    return true;
            }

            if(index == CYCLE_BASES.length - 1)
                index = 0;
            else
                ++index;
        }
    }

    private static int getIntegerAttribute(final SAMRecord record, final String tag, int defaultValue)
    {
        Object value = record.getAttribute(tag);
        return value != null ? (Integer)value : defaultValue;
    }
}
