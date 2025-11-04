package com.hartwig.hmftools.common.sequencing;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.PHRED_OFFSET;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;

import htsjdk.samtools.SAMRecord;

public final class UltimaBamUtils
{
    public static final byte ULTIMA_MAX_QUAL = 35;
    public static final int ULTIMA_MAX_HP_LEN = 20;

    public static final byte ULTIMA_INVALID_QUAL = -1;

    public static final String ULTIMA_TP_TAG = "tp";
    public static final String ULTIMA_T0_TAG = "t0";

    // PPM-SEQ tags
    private static final String PPM_STRAND_ST = "st";
    private static final String PPM_STRAND_ET = "et";
    private static final String PPM_TYPE_MIXED = "MIXED";

    public static final String ULT_QUAL_TAG = "UQ";
    public static final String ULT_QUAL_TAG_DELIM = "=";
    public static final String ULT_QUAL_TAG_INDEX_DELIM = ",";

    private static final int PPM_STRAND_MIN_SUM = 4;
    private static final int PPM_STRAND_MAX_SUM = 8;
    private static final double PPM_STRAND_BALANCED_LOW = 0.27;
    private static final double PPM_STRAND_BALANCED_HIGH = 1 - PPM_STRAND_BALANCED_LOW;

    // equivalent to adding their logs, ie P(combined) = 10^(-qual/10) + 10^(-qual/10), combined qual = -10 * log10(P(combined))
    public static final byte HALF_PHRED_SCORE_SCALING = 3;

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

    public static int extractLowQualCount(final SAMRecord record)
    {
        String qualTag = record.getStringAttribute(ULT_QUAL_TAG);

        if(qualTag == null)
            return 0;

        String[] qualItems = qualTag.split(ULT_QUAL_TAG_DELIM, 2);
        return Integer.parseInt(qualItems[0]);
    }

    public static List<Integer> extractLowQualIndices(final SAMRecord record)
    {
        String qualTag = record.getStringAttribute(ULT_QUAL_TAG);

        if(qualTag == null)
            return Collections.emptyList();

        String[] qualItems = qualTag.split(ULT_QUAL_TAG_DELIM, 2);

        int totalLowQual = Integer.parseInt(qualItems[0]);
        String[] qualIndexItems = qualItems[1].split(ULT_QUAL_TAG_INDEX_DELIM);

        List<Integer> lowQualValues = Lists.newArrayListWithCapacity(totalLowQual);

        for(String qualItem : qualIndexItems)
        {
            if(qualItem.contains("-"))
            {
                String[] startEnd = qualItem.split("-", 2);
                int indexStart = Integer.parseInt(startEnd[0]);
                int indexEnd = Integer.parseInt(startEnd[1]);

                for(int i = indexStart; i <= indexEnd; ++i)
                {
                    lowQualValues.add(i);
                }
            }
            else
            {
                lowQualValues.add(Integer.parseInt(qualItem));
            }
        }

        return lowQualValues;
    }

    public static ConsensusType deriveConsensusType(final SAMRecord record)
    {
        // st / et values: : MIXED, PLUS and MINUS
        String st = record.getStringAttribute(PPM_STRAND_ST);
        String et = record.getStringAttribute(PPM_STRAND_ET);

        if(st == null && et == null) // TEMP support for old tags
            return deriveConsensusTypeOld(record);

        if(st != null && et != null && st.equals(PPM_TYPE_MIXED) && et.equals(PPM_TYPE_MIXED))
            return ConsensusType.DUAL;

        return ConsensusType.NONE;
    }

    // old PPM-SEQ tags
    private static final String PPM_STRAND_AS = "as";
    private static final String PPM_STRAND_TS = "ts";
    private static final String PPM_STRAND_AE = "ae";
    private static final String PPM_STRAND_TE = "te";

    public static ConsensusType deriveConsensusTypeOld(final SAMRecord record)
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

    private static int getIntegerAttribute(final SAMRecord record, final String tag, int defaultValue)
    {
        Object value = record.getAttribute(tag);
        return value != null ? (Integer)value : defaultValue;
    }
}
