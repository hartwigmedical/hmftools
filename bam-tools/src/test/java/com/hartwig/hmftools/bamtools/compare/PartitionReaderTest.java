package com.hartwig.hmftools.bamtools.compare;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class PartitionReaderTest
{
    @Test
    public void testCompareReads()
    {
        final CompareConfig config = new CompareConfig();

        final SAMRecord read1 = new SAMRecord(null);
        read1.setReadName("read1");
        read1.setFirstOfPairFlag(true);
        read1.setReferenceName("chr1");
        read1.setAlignmentStart(2000);
        read1.setCigarString("151M");
        read1.setReadNegativeStrandFlag(false);
        read1.setMappingQuality(20);

        List<String> diffs = PartitionReader.compareReads(read1, read1, config);

        Assert.assertTrue(diffs.isEmpty());

        final SAMRecord read2 = new SAMRecord(null);
        read2.setReadName("read1");
        read2.setFirstOfPairFlag(true);
        read2.setReferenceName("chr1");
        read2.setAlignmentStart(2000);
        read2.setCigarString("151M");
        read2.setReadNegativeStrandFlag(false);
        read2.setMappingQuality(20);

        diffs = PartitionReader.compareReads(read1, read2, config);
        Assert.assertTrue(diffs.isEmpty());

        // different map qual
        read2.setMappingQuality(30);
        diffs = PartitionReader.compareReads(read1, read2, config);
        Assert.assertEquals(1, diffs.size());
        Assert.assertEquals("mapQuality(20/30)", diffs.get(0));

        // different read string
        read2.setMappingQuality(read1.getMappingQuality());
        read1.setReadString("ATCG");
        read2.setReadString("ATCC");
        diffs = PartitionReader.compareReads(read1, read2, config);
        Assert.assertEquals(1, diffs.size());
        Assert.assertEquals("bases(ATCG/ATCC)", diffs.get(0));

        // is negative strand
        read2.setReadNegativeStrandFlag(true);
        read1.setReadString("ATCG");
        read2.setReadString("CGAT");
        diffs = PartitionReader.compareReads(read1, read2, config);
        Assert.assertEquals(1, diffs.size());
        Assert.assertEquals("negStrand(false/true)", diffs.get(0));
    }
}