package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import java.util.List;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class PartitionReaderTest
{
    @Test
    public void testCompareReads()
    {
        CompareConfig config = new CompareConfig();

        String readId = "read1";
        String readBases = "ACGT";
        String readCigar = "4M";

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 2000, readBases, readCigar, CHR_1, 2200,
                false, false, null);
        read1.setMappingQuality(20);

        List<String> diffs = PartitionReader.compareReads(read1, read1, config);

        Assert.assertTrue(diffs.isEmpty());

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 2000, readBases, readCigar, CHR_1, 2200,
                false, false, null);
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