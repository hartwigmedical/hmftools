package com.hartwig.hmftools.common.sequencing;

import htsjdk.samtools.SAMRecord;

public final class UltimaBamUtils
{
    public static final byte ULTIMA_MAX_QUAL = 40;

    private static final String PPM_STRAND_AS = "as";
    private static final String PPM_STRAND_TS = "ts";
    private static final String PPM_STRAND_AE = "ae";
    private static final String PPM_STRAND_TE = "te";

    private static final int PPM_STRAND_MIN_SUM = 4;
    private static final int PPM_STRAND_MAX_SUM = 8;
    private static final double PPM_STRAND_BALANCED_LOW = 0.27;
    private static final double PPM_STRAND_BALANCED_HIGH = 1 - PPM_STRAND_BALANCED_LOW;

    private static final String TP_TAG = "tp";
    private static final String T0_TAG = "t0";

    public static byte[] extractTpValues(final SAMRecord record)
    {
        byte[] results = new byte[record.getBaseQualities().length];

        byte[] rawValues = record.getByteArrayAttribute(TP_TAG);

        return results;
    }

    public static short[] extractT0Values(final SAMRecord record)
    {
        String strValue = record.getStringAttribute(T0_TAG);
        short[] results = new short[record.getBaseQualities().length];
        return results;
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

    private static int getIntegerAttribute(final SAMRecord record, final String tag, int defaultValue)
    {
        Object value = record.getAttribute(tag);
        return value != null ? (Integer)value : defaultValue;
    }
}
