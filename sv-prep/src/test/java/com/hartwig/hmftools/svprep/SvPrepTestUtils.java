package com.hartwig.hmftools.svprep;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public final class SvPrepTestUtils
{
    public static final String CHR_1 = "1";
    public static final int DEFAULT_MAP_QUAL = 60;
    public static final int DEFAULT_BASE_QUAL = 37;

    public static int buildFlags(boolean firstInPair, boolean reversed, boolean duplicate, boolean supplementary, boolean secondary)
    {
        int flags = 0;

        flags = setReadFlag(flags, SAMFlag.READ_PAIRED);
        flags = setReadFlag(flags, SAMFlag.PROPER_PAIR);

        if(reversed)
            flags = setReadFlag(flags, SAMFlag.READ_REVERSE_STRAND);

        if(firstInPair)
            flags = setReadFlag(flags, SAMFlag.FIRST_OF_PAIR);
        else
            flags = setReadFlag(flags, SAMFlag.SECOND_OF_PAIR);

        if(secondary)
            flags = setReadFlag(flags, SAMFlag.SECONDARY_ALIGNMENT);

        if(supplementary)
            flags = setReadFlag(flags, SAMFlag.SUPPLEMENTARY_ALIGNMENT);

        return flags;
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chromosome, int readStart, final String readBases, final String cigar)
    {
        return createSamRecord(
                readId, chromosome, readStart, readBases, cigar,
                buildFlags(true, false, false, false, false),
                DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL);
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chromosome, int readStart, final String readBases, final String cigar, int flags,
            int mapQual, int baseQual)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);
        SAMRecord record = recordBuilder.addFrag(
                readId, 1, readStart, false, false, cigar, readBases, mapQual, false);

        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = (byte)baseQual;

        record.setBaseQualities(qualities);
        record.setReferenceName(chromosome);

        record.setFlags(flags);
        return record;
    }

    public static int setReadFlag(int flags, final SAMFlag flag)
    {
        flags |= flag.intValue();
        return flags;
    }
}
