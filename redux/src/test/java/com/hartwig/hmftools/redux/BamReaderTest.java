package com.hartwig.hmftools.redux;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class BamReaderTest
{
    private static SAMRecord mappedRead(final Cigar cigar, final int seqLength)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read1");
        record.setReadUnmappedFlag(false);
        record.setCigar(cigar);

        byte[] bases = new byte[seqLength];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        return record;
    }

    private static Cigar cigar(final CigarElement... elements)
    {
        return new Cigar(Arrays.asList(elements));
    }

    @Test
    public void testValidReadNotMalformed()
    {
        assertFalse(BamReader.isMalformed(mappedRead(cigar(new CigarElement(10, CigarOperator.M)), 10)));
    }

    @Test
    public void testSeqShorterThanCigarIsMalformed()
    {
        assertTrue(BamReader.isMalformed(mappedRead(cigar(new CigarElement(10, CigarOperator.M)), 5)));
    }

    @Test
    public void testZeroLengthCigarElementIsMalformed()
    {
        assertTrue(BamReader.isMalformed(mappedRead(
                cigar(new CigarElement(0, CigarOperator.M), new CigarElement(10, CigarOperator.M)), 10)));
    }

    @Test
    public void testEmptyCigarNotMalformed()
    {
        assertFalse(BamReader.isMalformed(mappedRead(new Cigar(), 10)));
    }

    @Test
    public void testUnmappedReadNotMalformed()
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadUnmappedFlag(true);
        assertFalse(BamReader.isMalformed(record));
    }
}
