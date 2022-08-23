package com.hartwig.hmftools.common.samtools;

import static com.hartwig.hmftools.common.samtools.CigarUtils.calcCigarLength;

import static org.junit.Assert.assertEquals;
import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

public class CigarUtilsTest
{
    @Test
    public void testCigarFromStr()
    {
        Cigar cigar = CigarUtils.cigarFromStr("120S35M1099N11M5H");

        assertEquals(5, cigar.numCigarElements());
        assertEquals(120, cigar.getCigarElement(0).getLength());
        assertEquals(CigarOperator.S, cigar.getCigarElement(0).getOperator());
        assertEquals(35, cigar.getCigarElement(1).getLength());
        assertEquals(CigarOperator.M, cigar.getCigarElement(1).getOperator());
        assertEquals(1099, cigar.getCigarElement(2).getLength());
        assertEquals(CigarOperator.N, cigar.getCigarElement(2).getOperator());
        assertEquals(11, cigar.getCigarElement(3).getLength());
        assertEquals(CigarOperator.M, cigar.getCigarElement(3).getOperator());
        assertEquals(5, cigar.getCigarElement(4).getLength());
        assertEquals(CigarOperator.H, cigar.getCigarElement(4).getOperator());
    }

    @Test
    public void testCalcCigarLength()
    {
        TestCase.assertEquals(100, calcCigarLength("100M"));
        TestCase.assertEquals(100, calcCigarLength("10M10D80M"));
        TestCase.assertEquals(100, calcCigarLength("10M1I1000N10D2I80M"));
        TestCase.assertEquals(10, calcCigarLength("9M10I1D"));

        // error handling
        TestCase.assertEquals(0, calcCigarLength("12345"));
        TestCase.assertEquals(0, calcCigarLength("G10Z2QR100"));
    }

    @Test
    public void testSoftClip()
    {
        // left soft clip
        SAMRecord record = new SAMRecord(null);
        record.setCigarString("10S85M");
        TestCase.assertEquals(10, CigarUtils.leftSoftClip(record));
        TestCase.assertEquals(0, CigarUtils.rightSoftClip(record));

        // right soft clip
        record = new SAMRecord(null);
        record.setCigarString("85M10S");
        TestCase.assertEquals(0, CigarUtils.leftSoftClip(record));
        TestCase.assertEquals(10, CigarUtils.rightSoftClip(record));

        // hard clip should not count
        record = new SAMRecord(null);
        record.setCigarString("10H85M7H");
        TestCase.assertEquals(0, CigarUtils.leftSoftClip(record));
        TestCase.assertEquals(0, CigarUtils.rightSoftClip(record));
    }
}
