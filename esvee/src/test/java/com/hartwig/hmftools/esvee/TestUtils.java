package com.hartwig.hmftools.esvee;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.read.Read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class TestUtils
{
    public static final String REF_BASES = generateRandomBases(100);
    public static final String TEST_READ_ID = "READ_01";
    public static final String TEST_CIGAR  = "100M";

    public static Read createSAMRecord(final String sequence)
    {
        return createSAMRecord(sequence, 1);
    }

    public static Read createSAMRecord(final String sequence, final int position)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReferenceName(CHR_1);
        record.setReadName(String.valueOf(System.nanoTime()));
        record.setCigarString(sequence.length() + "M");
        record.setAlignmentStart(position);
        record.setReadString(sequence);
        record.setBaseQualityString("F".repeat(sequence.length()));
        record.setReadPairedFlag(true);
        record.setMateUnmappedFlag(true);
        return new Read(record);
    }

    public static Read createSamRecord(final String readId, int readStart, final String readBases, final String cigar)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, readStart + 1000,
                false, false, null);
        return new Read(record);
    }

    public static String makeCigarString(final String readBases, int scLeft, int scRight)
    {
        StringBuilder sb = new StringBuilder();

        if(scLeft > 0)
            sb.append(format("%dS", scLeft));

        sb.append(format("%dM", readBases.length() - scLeft - scRight));

        if(scRight > 0)
            sb.append(format("%dS", scRight));

        return sb.toString();
    }
}