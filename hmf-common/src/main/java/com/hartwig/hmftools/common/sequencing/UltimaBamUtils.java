package com.hartwig.hmftools.common.sequencing;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.PHRED_OFFSET;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;

import htsjdk.samtools.SAMRecord;

public final class UltimaBamUtils
{
    public static final byte ULTIMA_MAX_QUAL_TP = 40;
    public static final byte TP_0_BOOST = 5;
    public static final byte ULTIMA_MAX_QUAL_T0 = 40;
    public static final byte BOOSTED_QUAL = 35;

    public static final byte ULTIMA_INVALID_QUAL = -1;

    public static final String TP_TAG = "tp";
    public static final String T0_TAG = "t0";

    private static final String PPM_STRAND_AS = "as";
    private static final String PPM_STRAND_TS = "ts";
    private static final String PPM_STRAND_AE = "ae";
    private static final String PPM_STRAND_TE = "te";

    private static final int PPM_STRAND_MIN_SUM = 4;
    private static final int PPM_STRAND_MAX_SUM = 8;
    private static final double PPM_STRAND_BALANCED_LOW = 0.27;
    private static final double PPM_STRAND_BALANCED_HIGH = 1 - PPM_STRAND_BALANCED_LOW;

    public static final byte[] CYCLE_BASES = new byte[] { (byte) 'T', (byte) 'G', (byte) 'C', (byte) 'A' };

    public static byte[] extractTpValues(final SAMRecord record)
    {
        return record.getByteArrayAttribute(TP_TAG);
    }

    public static byte[] extractT0Values(final SAMRecord record)
    {
        byte[] t0Values = record.getStringAttribute(T0_TAG).getBytes();
        for(int i = 0; i < t0Values.length; ++i)
            t0Values[i] -= PHRED_OFFSET;
        return t0Values;
    }

    public static UltimaConsensusType extractConsensusType(final SAMRecord record)
    {
        //  example: as:i:3  ts:i:4  ae:i:3  te:i:3
        int as = getIntegerAttribute(record, PPM_STRAND_AS, 0);
        int ts = getIntegerAttribute(record, PPM_STRAND_TS, 0);
        int ae = getIntegerAttribute(record, PPM_STRAND_AE, 0);
        int te = getIntegerAttribute(record, PPM_STRAND_TE, 0);

        double sumStart = as + ts;
        double sumEnd = ae + te;

        if(sumStart < PPM_STRAND_MIN_SUM || sumStart > PPM_STRAND_MAX_SUM || sumEnd < PPM_STRAND_MIN_SUM || sumEnd > PPM_STRAND_MAX_SUM)
            return UltimaConsensusType.STANDARD;

        double ratioStart = as / sumStart;
        double ratioEnd = ae / sumEnd;

        if(ratioStart < PPM_STRAND_BALANCED_LOW || ratioStart > PPM_STRAND_BALANCED_HIGH)
            return UltimaConsensusType.STANDARD;

        if(ratioEnd < PPM_STRAND_BALANCED_LOW || ratioEnd > PPM_STRAND_BALANCED_HIGH)
            return UltimaConsensusType.STANDARD;

        return UltimaConsensusType.BALANCED;
    }

    public static byte calcTpBaseQual(final SAMRecord record, int indexStart, int indexEnd, int tpSearchValue)
    {
        if(indexStart < 0)
        {
            return ULTIMA_INVALID_QUAL;
        }

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
        if(qualValue1 == BOOSTED_QUAL)
            qualValue1 = ULTIMA_MAX_QUAL_TP;
        if(qualValue2 == BOOSTED_QUAL)
            qualValue2 = ULTIMA_MAX_QUAL_TP;

        if(qualValue1 < 0)
        {
            // if bases aren't found, use middle base(s) as: min(40, Q + 6 * abs(required_tp – T))
            int middleIndex = indexStart + (indexEnd - indexStart) / 2;
            final byte[] baseQualities = record.getBaseQualities();
            if(middleIndex < 0 || middleIndex >= baseQualities.length)
            {
                return ULTIMA_INVALID_QUAL;
            }

            qualValue1 = baseQualities[middleIndex];
            byte tpValue = tpValues[middleIndex];
            if(tpValue == (byte) 0)
                return ULTIMA_MAX_QUAL_TP + TP_0_BOOST;

            if(qualValue1 == BOOSTED_QUAL)
                qualValue1 = ULTIMA_MAX_QUAL_TP;

            return (byte)min(ULTIMA_MAX_QUAL_TP, qualValue1 + 6 * abs(tpSearchValue - tpValue));
        }

        if(qualValue2 < 0)
            return (byte)qualValue1;

        // equivalent to adding their logs, ie P(combined) = 10^(-qual/10) + 10^(-qual/10), combined qual = -10 * log10(P(combined))
        return (byte)(qualValue1 - 3);
    }


    public static boolean isBaseInCycle(final byte startBase, final byte endBase, final byte testBase)
    {
        if(startBase == endBase || startBase == testBase)
            return false;

        if(baseIndex(startBase) == -1 || baseIndex(endBase) == -1)
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

    private static int getIntegerAttribute(final SAMRecord record, final String tag, int defaultValue)
    {
        Object value = record.getAttribute(tag);
        return value != null ? (Integer)value : defaultValue;
    }

    public static byte safeQualLookup(final byte[] qualityArray, int index)
    {
        if(index < 0 || index >= qualityArray.length)
            return ULTIMA_INVALID_QUAL;

        byte value = qualityArray[index];
        if(value == BOOSTED_QUAL)
            return ULTIMA_MAX_QUAL_T0;

        return value;
    }
}
