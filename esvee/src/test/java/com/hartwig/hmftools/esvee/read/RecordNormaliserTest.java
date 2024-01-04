package com.hartwig.hmftools.esvee.read;

import java.util.Arrays;

import com.hartwig.hmftools.common.samtools.CigarUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class RecordNormaliserTest
{
    private Read record(final String cigarString)
    {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("1", 100_000));
        final SAMRecord record = new SAMRecord(header);

        final Cigar cigar = CigarUtils.cigarFromStr(cigarString);
        final int length = cigar.getReadLength();
        final String bases = "A".repeat(length);
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 37);

        record.setReadString(bases);
        record.setBaseQualities(quals);
        record.setReferenceName("1");
        record.setAlignmentStart(1000);
        record.setCigar(cigar);

        return new Read(record);
    }

    /* CHASHA FIXME
    @Test
    public void handlesMalformedCigars()
    {
        final RecordNormaliser normaliser = new RecordNormaliser(new MockRefGenome(), TestUtils.config());

        final Record record = record("15M10I100S");
        final Record result = normaliser.normalise(record);

        assertTrue(record).isNotSameAs(result);
        assertTrue(result).isNotNull();
        assertTrue(result.getCigar().toString()).isEqualTo("15M110S");
    }

     */
}