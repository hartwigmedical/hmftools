package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import com.hartwig.hmftools.svprep.reads.ReadFilterConfig;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public final class SvPrepTestUtils
{
    public static final String CHR_1 = "1";
    public static final String CHR_2 = "2";
    public static final String CHR_3 = "3";
    public static final String CHR_4 = "4";

    public static final int DEFAULT_MAP_QUAL = 60;
    public static final int DEFAULT_BASE_QUAL = 37;

    public static final ReadFilterConfig READ_FILTERS = ReadFilterConfig.from(null);
    public static final HotspotCache HOTSPOT_CACHE = new HotspotCache(null);
    public static final BlacklistLocations BLACKLIST_LOCATIONS = new BlacklistLocations(null);

    public static String readIdStr(int readId) { return format("READ_%02d", readId); }

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
            final String readId, final String chromosome, int readStart, final String mateChr, int mateStart,
            boolean firstInPair, boolean isSupp, final String suppData)
    {
        SAMRecord record = createSamRecord(
                readId, chromosome, readStart, "", "100M",
                buildFlags(firstInPair, false, false, isSupp, false),
                DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL);

        record.setMateReferenceName(mateChr);
        record.setMateAlignmentStart(mateStart);

        if(!suppData.isEmpty())
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, suppData);

        return record;
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
