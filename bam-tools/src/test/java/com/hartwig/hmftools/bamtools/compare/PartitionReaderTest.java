package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import java.util.List;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class PartitionReaderTest
{
    private static final String TEST_READ_ID = "read1";
    private static final String TEST_READ_BASES = "ACGT";
    private static final String TEST_CIGAR = "4M";

    @Test
    public void testCompareReads()
    {
        CompareConfig config = new CompareConfig();

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_1, 2000, TEST_READ_BASES, TEST_CIGAR, CHR_1, 2200,
                false, false, null);
        read1.setMappingQuality(20);

        List<String> diffs = CompareUtils.compareReads(read1, read1, config);

        Assert.assertTrue(diffs.isEmpty());

        SAMRecord read2 = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_1, 2000, TEST_READ_BASES, TEST_CIGAR, CHR_1, 2200,
                false, false, null);
        read2.setMappingQuality(20);

        diffs = CompareUtils.compareReads(read1, read2, config);
        Assert.assertTrue(diffs.isEmpty());

        // different map qual
        read2.setMappingQuality(30);
        diffs = CompareUtils.compareReads(read1, read2, config);
        Assert.assertEquals(1, diffs.size());
        Assert.assertEquals("mapQuality(20/30)", diffs.get(0));

        // different read string
        read2.setMappingQuality(read1.getMappingQuality());

        read1.setReadString("ATCG");
        read2.setReadString("ATCC");
        diffs = CompareUtils.compareReads(read1, read2, config);
        Assert.assertEquals(1, diffs.size());
        Assert.assertEquals("bases(ATCG/ATCC)", diffs.get(0));

        // is negative strand
        read2.setReadNegativeStrandFlag(true);
        read1.setReadString("ATCG");
        read2.setReadString("CGAT");
        diffs = CompareUtils.compareReads(read1, read2, config);
        Assert.assertEquals(1, diffs.size());
        Assert.assertEquals("negStrand(false/true)", diffs.get(0));
    }

    @Test
    public void testReduxReadAdjustments()
    {
        CompareConfig config = new CompareConfig(
                true, true, true,
                 true, true, true);

        // reads differ on duplicate status
        int readPosition = 1000;
        int matePosition = 2000;

        SAMRecord read1 = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_1, readPosition, TEST_READ_BASES, TEST_CIGAR, CHR_2, matePosition,
                false, false, null, true, TEST_CIGAR);

        SAMRecord read2 = SamRecordTestUtils.cloneSamRecord(read1, TEST_READ_ID);

        read2.setDuplicateReadFlag(true);

        List<String> diffs = CompareUtils.compareReads(read1, read2, config);
        Assert.assertTrue(diffs.isEmpty());

        // read is unmapped
        read2 = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, NO_CHROMOSOME_NAME, matePosition, TEST_READ_BASES, NO_CIGAR, CHR_2, matePosition,
                false, false, null);

        read2.setReadUnmappedFlag(true);
        read2.setInferredInsertSize(-1);

        // no mate cigar set

        diffs = CompareUtils.compareReads(read1, read2, config);
        Assert.assertTrue(diffs.isEmpty());


    }
}