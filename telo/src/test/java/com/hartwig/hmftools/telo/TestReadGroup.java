package com.hartwig.hmftools.telo;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import junit.framework.TestCase;

import java.util.List;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.junit.Test;

public class TestReadGroup
{
    @Test
    public void testSuppAlignmentPositions()
    {
        String saString = "16,59998,+,64S38M49S,55,0;";

        List<ChrBaseRegion> brs = ReadGroup.suppAlignmentPositions(saString);

        TestCase.assertEquals(brs.size(), 1);
        TestCase.assertEquals(brs.get(0).Chromosome, "16");
        TestCase.assertEquals(brs.get(0).start(), 59998);

        saString = "16,59998,+,64S38M49S,55,0;15,42243201,+,108S35M8S,9,0;,";

        brs = ReadGroup.suppAlignmentPositions(saString);

        TestCase.assertEquals(brs.size(), 2);
        TestCase.assertEquals(brs.get(0).Chromosome, "16");
        TestCase.assertEquals(brs.get(0).start(), 59998);

        TestCase.assertEquals(brs.get(1).Chromosome, "15");
        TestCase.assertEquals(brs.get(1).start(), 42243201);
    }

    @Test
    public void testReadGroupIsComplete()
    {
        ReadGroup rg = new ReadGroup("ReadGroup");

        SAMRecord record = new SAMRecord(null);
        record.setReadName(rg.getName());
        record.setReferenceName("1");
        record.setAlignmentStart(100);
        record.setMateReferenceName("2");
        record.setMateAlignmentStart(200);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        rg.Reads.add(record);

        TestCase.assertFalse(rg.isComplete());

        record = new SAMRecord(null);
        record.setReadName(rg.getName());
        record.setReferenceName("2");
        record.setAlignmentStart(200);
        record.setMateReferenceName("1");
        record.setMateAlignmentStart(100);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(false);
        rg.Reads.add(record);

        TestCase.assertTrue(rg.isComplete());
        TestCase.assertTrue(rg.invariant());
    }

    @Test
    public void testReadGroupIsCompleteWithSupplementary()
    {
        ReadGroup rg = new ReadGroup("ReadGroup");

        SAMRecord record = new SAMRecord(null);
        record.setReadName(rg.getName());
        record.setReferenceName("1");
        record.setAlignmentStart(100);
        record.setMateReferenceName("2");
        record.setMateAlignmentStart(200);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        rg.Reads.add(record);

        TestCase.assertFalse(rg.isComplete());

        record = new SAMRecord(null);
        record.setReadName(rg.getName());
        record.setReferenceName("2");
        record.setAlignmentStart(200);
        record.setMateReferenceName("1");
        record.setMateAlignmentStart(100);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(false);

        record.setAttribute(SAMTag.SA.name(), "3,300,+,64S38M49S,55,0;4,400,+,64S38M49S,55,0;,");

        rg.Reads.add(record);

        TestCase.assertFalse(rg.isComplete());

        // add supplementary
        record = new SAMRecord(null);
        record.setReadName(rg.getName());
        record.setReferenceName("3");
        record.setAlignmentStart(300);
        record.setMateReferenceName("1");
        record.setMateAlignmentStart(100);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(false);
        record.setSupplementaryAlignmentFlag(true);

        record.setAttribute(SAMTag.SA.name(), "2,200,+,64S38M49S,55,0;");

        rg.SupplementaryReads.add(record);

        // we still missing one supplementary read
        TestCase.assertFalse(rg.isComplete());

        record = new SAMRecord(null);
        record.setReadName(rg.getName());
        record.setReferenceName("4");
        record.setAlignmentStart(400);
        record.setMateReferenceName("1");
        record.setMateAlignmentStart(100);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(false);
        record.setSupplementaryAlignmentFlag(true);

        record.setAttribute(SAMTag.SA.name(), "2,200,+,64S38M49S,55,0;");

        rg.SupplementaryReads.add(record);

        TestCase.assertTrue(rg.isComplete());
        TestCase.assertTrue(rg.invariant());
    }
}
