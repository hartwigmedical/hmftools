package com.hartwig.hmftools.bamtools.compare;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

public class ReadKeyTest
{
    // test that we can create keys that differentiate supplementary reads
    @Test
    public void testCalcSupplementaryIndex()
    {
        final SAMRecord read = new SAMRecord(null);
        read.setReadName("read1");
        read.setFirstOfPairFlag(true);
        read.setReferenceName("chr1");
        read.setAlignmentStart(2000);
        read.setCigarString("151M");
        read.setReadNegativeStrandFlag(false);
        read.setMappingQuality(20);

        // this is not a supplementary, should have 0 as supplementary index
        ReadKey readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(0, readKey.SupplementaryIndex);

        // add supplementary attribute
        read.setAttribute(SAMTag.SA.name(), "chr1,2000,+,151M,0,1;");

        // still since it is not supplementary, should have 0 as supplementary index
        readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(0, readKey.SupplementaryIndex);

        // make it a supplementary
        read.setSupplementaryAlignmentFlag(true);

        // The SA attribute shows only the non supplementary read plus this read, so the supplementary index of this
        // read should be 1
        readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(1, readKey.SupplementaryIndex);

        // set SA to show 2 supplementary reads, this time the other supplementary read is mapped to chr1:3000
        // which is > this read's genomic location, as a result the supplementary index is still 1
        read.setAttribute(SAMTag.SA.name(), "chr1,1000,+,151M,0,1;chr1,3000,+,151M,0,1;");

        readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(1, readKey.SupplementaryIndex);

        // set SA to show 2 supplementary reads, this time the other supplementary read is mapped to chr1:1000
        // which is < this read's genomic location, as a result the supplementary index will be 2
        read.setAttribute(SAMTag.SA.name(), "chr1,1000,+,151M,0,1;chr1,1000,+,151M,0,1;");

        readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(2, readKey.SupplementaryIndex);

        // set SA to show 3 supplementary reads, both are located genomically lower than this read
        // As a result the supplementary index will be 3
        read.setAttribute(SAMTag.SA.name(), "chr1,1000,+,151M,0,1;chr1,1000,+,151M,0,1;chr1,1500,+,151M,0,1;");

        readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(3, readKey.SupplementaryIndex);

        // set SA to show 3 supplementary reads, both are located genomically higher than this read
        // As a result the supplementary index will be 1
        read.setAttribute(SAMTag.SA.name(), "chr1,3000,+,151M,0,1;chr1,4000,+,151M,0,1;chr1,5000,+,151M,0,1;");

        readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(1, readKey.SupplementaryIndex);

        // set SA to show 3 supplementary reads, first one is located genomically higher than this read, 2nd one is
        // located genomically lower.
        // As a result the supplementary index will be 2
        read.setAttribute(SAMTag.SA.name(), "chr1,1000,+,151M,0,1;chr1,3000,+,151M,0,1;chr1,1000,+,151M,0,1;");

        readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(2, readKey.SupplementaryIndex);

        // set SA to show 3 supplementary reads, first one is located genomically lower than this read, 2nd one is
        // located genomically higher.
        // As a result the supplementary index will be 2
        read.setAttribute(SAMTag.SA.name(), "chr1,1000,+,151M,0,1;chr1,1000,+,151M,0,1;chr1,3000,+,151M,0,1;");

        readKey = ReadKey.from(read);
        Assert.assertEquals("read1", readKey.ReadName);
        Assert.assertTrue(readKey.FirstOfPair);
        Assert.assertEquals(2, readKey.SupplementaryIndex);
    }
}