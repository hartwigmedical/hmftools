package com.hartwig.hmftools.common.bam;

import static com.hartwig.hmftools.common.bam.CigarUtils.calcCigarAlignedLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.getPositionFromReadIndex;
import static com.hartwig.hmftools.common.bam.CigarUtils.getReadIndexFromPosition;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.INVALID_READ_INDEX;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
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
        TestCase.assertEquals(100, calcCigarAlignedLength("100M"));
        TestCase.assertEquals(100, calcCigarAlignedLength("10M10D80M"));
        TestCase.assertEquals(1100, calcCigarAlignedLength("10M1I1000N10D2I80M"));
        TestCase.assertEquals(10, calcCigarAlignedLength("9M10I1D"));

        // error handling
        TestCase.assertEquals(0, calcCigarAlignedLength("12345"));
        TestCase.assertEquals(0, calcCigarAlignedLength("G10Z2QR100"));
    }

    @Test
    public void testSoftClip()
    {
        // left soft clip
        SAMRecord record = new SAMRecord(null);
        record.setCigarString("10S85M");
        TestCase.assertEquals(10, CigarUtils.leftSoftClipLength(record));
        TestCase.assertEquals(0, CigarUtils.rightSoftClipLength(record));

        // right soft clip
        record = new SAMRecord(null);
        record.setCigarString("85M10S");
        TestCase.assertEquals(0, CigarUtils.leftSoftClipLength(record));
        TestCase.assertEquals(10, CigarUtils.rightSoftClipLength(record));

        // hard clip should not count
        record = new SAMRecord(null);
        record.setCigarString("10H85M7H");
        TestCase.assertEquals(0, CigarUtils.leftSoftClipLength(record));
        TestCase.assertEquals(0, CigarUtils.rightSoftClipLength(record));
    }

    @Test
    public void testReadIndexFromPosition()
    {
        List<CigarElement> cigarElements = Lists.newArrayList();

        cigarElements.add(new CigarElement(20, CigarOperator.M));

        int readIndex = getReadIndexFromPosition(100, cigarElements, 100);
        assertEquals(0, readIndex);

        readIndex = getReadIndexFromPosition(100, cigarElements, 119);
        assertEquals(19, readIndex);

        // outisde the bounds
        readIndex = getReadIndexFromPosition(100, cigarElements, 120);
        assertEquals(INVALID_READ_INDEX, readIndex);

        readIndex = getReadIndexFromPosition(100, cigarElements, 99);
        assertEquals(INVALID_READ_INDEX, readIndex);

        // within a delete
        cigarElements.clear();
        cigarElements.add(new CigarElement(20, CigarOperator.M));
        cigarElements.add(new CigarElement(10, CigarOperator.D));
        cigarElements.add(new CigarElement(20, CigarOperator.M));

        readIndex = getReadIndexFromPosition(100, cigarElements, 125, false, false);
        assertEquals(INVALID_READ_INDEX, readIndex);

        readIndex = getReadIndexFromPosition(100, cigarElements, 125, true, false);
        assertEquals(19, readIndex);

        // within soft-clips
        cigarElements.clear();
        cigarElements.add(new CigarElement(10, CigarOperator.S));
        cigarElements.add(new CigarElement(20, CigarOperator.M));
        cigarElements.add(new CigarElement(10, CigarOperator.S));

        readIndex = getReadIndexFromPosition(100, cigarElements, 95, false, false);
        assertEquals(INVALID_READ_INDEX, readIndex);

        readIndex = getReadIndexFromPosition(100, cigarElements, 95, false, true);
        assertEquals(5, readIndex);

        readIndex = getReadIndexFromPosition(100, cigarElements, 125, false, false);
        assertEquals(INVALID_READ_INDEX, readIndex);

        readIndex = getReadIndexFromPosition(100, cigarElements, 125, false, true);
        assertEquals(35, readIndex);
    }

    @Test
    public void testPositionFromReadIndex()
    {
        List<CigarElement> cigarElements = Lists.newArrayList();

        cigarElements.add(new CigarElement(20, CigarOperator.M));

        int readIndex = getPositionFromReadIndex(100, cigarElements, 0);
        assertEquals(100, readIndex);

        readIndex = getPositionFromReadIndex(100, cigarElements, 19);
        assertEquals(119, readIndex);

        readIndex = getPositionFromReadIndex(100, cigarElements, -1);
        assertEquals(NO_POSITION, readIndex);

        readIndex = getPositionFromReadIndex(100, cigarElements, 20);
        assertEquals(NO_POSITION, readIndex);

        // within a insert
        cigarElements.clear();
        cigarElements.add(new CigarElement(20, CigarOperator.M));
        cigarElements.add(new CigarElement(10, CigarOperator.I));
        cigarElements.add(new CigarElement(20, CigarOperator.M));

        readIndex = getPositionFromReadIndex(100, cigarElements, 25);
        assertEquals(NO_POSITION, readIndex);

        int[] readInfo = getPositionFromReadIndex(
                100, cigarElements, 25, true, false);

        assertEquals(119, readInfo[0]);
        assertEquals(-6, readInfo[1]);

        // within soft-clips
        cigarElements.clear();
        cigarElements.add(new CigarElement(10, CigarOperator.S));
        cigarElements.add(new CigarElement(20, CigarOperator.M));
        cigarElements.add(new CigarElement(10, CigarOperator.S));

        readIndex = getPositionFromReadIndex(100, cigarElements, 5);
        assertEquals(NO_POSITION, readIndex);

        readIndex = getPositionFromReadIndex(100, cigarElements, 35);
        assertEquals(NO_POSITION, readIndex);

        readIndex = getPositionFromReadIndex(100, cigarElements, 5, false, true)[0];
        assertEquals(95, readIndex);

        readIndex = getPositionFromReadIndex(100, cigarElements, 35, false, true)[0];
        assertEquals(125, readIndex);
    }
}
