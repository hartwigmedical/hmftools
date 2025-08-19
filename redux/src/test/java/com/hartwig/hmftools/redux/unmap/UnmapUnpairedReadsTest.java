package com.hartwig.hmftools.redux.unmap;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.ALIGNMENTS_DELIM;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.redux.TestUtils.checkTransformRead;
import static com.hartwig.hmftools.redux.unmap.UnmapReadsTest.READ_BASES;
import static com.hartwig.hmftools.redux.unmap.UnmapReadsTest.SOFT_CLIPPED_READ_CIGAR;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UnmapUnpairedReadsTest
{
    @Test
    public void testNoUnmapUnpaired()
    {
        // read not in unmap regions
        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, 100, READ_BASES, TEST_READ_CIGAR, false, false, null);

        assertFalse(checkTransformRead(read, CHR_1));
        assertFalse(read.getReadUnmappedFlag());

    }

    @Test
    public void testNoUnmapSoftClippedUnpaired()
    {
        // soft-clipped but not in an unmap region
        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, 100, READ_BASES, SOFT_CLIPPED_READ_CIGAR, false, false, null);

        assertFalse(checkTransformRead(read, CHR_1));
        assertFalse(read.getReadUnmappedFlag());
    }

    @Test
    public void testUnmapReadSoftClippedUnpaired()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, 550, READ_BASES, SOFT_CLIPPED_READ_CIGAR, false, false, null);

        assertTrue(checkTransformRead(read, CHR_1));
        assertTrue(read.getReadUnmappedFlag());
        assertEquals(NO_POSITION, read.getAlignmentStart());
        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
    }

    @Test
    public void testSupplementarySoftClippedUnpaired()
    {
        // supp removed since it is in an unmap region and has a long SC
        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, 100, READ_BASES, TEST_READ_CIGAR, false,
                false,
                new SupplementaryReadData(CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60));

        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(checkTransformRead(read, CHR_1));
        assertFalse(read.getReadUnmappedFlag());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapSupplementaryWhenPrimaryIsUnmappedSoftClippedUnpaired()
    {
        // Note that the suppReadData for a supplementary read refers to the associated primary read
        SupplementaryReadData suppReadData = new SupplementaryReadData(
                CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60);

        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_2, 100, READ_BASES, TEST_READ_CIGAR, false,
                true, suppReadData);

        // when supplementaries are unmapped not all properties are unset since it will be dropped immediately
        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(checkTransformRead(read, CHR_2));
        assertTrue(read.getReadUnmappedFlag());
    }

    @Test
    public void testUnpairedChimericNotUnmapped()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                TEST_READ_ID, CHR_1, 550, READ_BASES, TEST_READ_CIGAR, CHR_2, 600, false,
                false, null, true, TEST_READ_CIGAR);
        read.setReadPairedFlag(false);

        assertFalse(checkTransformRead(read, CHR_1));
        assertFalse(read.getReadUnmappedFlag());
    }

    @Test
    public void testUnmapReadHighDepthUnpaired()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(TEST_READ_ID, CHR_3, 600, READ_BASES, TEST_READ_CIGAR, false, false, null);

        assertTrue(checkTransformRead(read, CHR_3));
        assertTrue(read.getReadUnmappedFlag());
        assertEquals(NO_POSITION, read.getAlignmentStart());
        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
    }

    @Test
    public void testUnmapSupplementaryHighDepthUnpaired()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_3, 400, READ_BASES, TEST_READ_CIGAR, false,
                false,
                new SupplementaryReadData(CHR_3, 600, SUPP_POS_STRAND, TEST_READ_CIGAR, 60));

        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(checkTransformRead(read, CHR_3));
        assertFalse(read.getReadUnmappedFlag());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapMultipleSupplementaryAlignmentsUnpaired()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, 100, READ_BASES, TEST_READ_CIGAR, false,
                false, null);

        SupplementaryReadData suppData1 = new SupplementaryReadData(CHR_1, 1000, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
        SupplementaryReadData suppData2 = new SupplementaryReadData(CHR_2, 6000, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
        SupplementaryReadData suppData3 = new SupplementaryReadData(CHR_3, 600, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

        StringJoiner saJoiner = new StringJoiner(ALIGNMENTS_DELIM);
        saJoiner.add(suppData1.asSamTag());
        saJoiner.add(suppData2.asSamTag());
        saJoiner.add(suppData3.asSamTag());

        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, saJoiner.toString());

        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(checkTransformRead(read, CHR_1));
        assertFalse(read.getReadUnmappedFlag());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapSupplementaryWhenUnmappingAssociatedSupplementaryUnpaired()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_1, 100, READ_BASES, TEST_READ_CIGAR, false,
                true, null);

        SupplementaryReadData suppDataPrimary = new SupplementaryReadData(CHR_3, 100, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);
        SupplementaryReadData suppDataOther = new SupplementaryReadData(CHR_3, 500, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

        StringJoiner saJoiner = new StringJoiner(ALIGNMENTS_DELIM);
        saJoiner.add(suppDataPrimary.asSamTag());
        saJoiner.add(suppDataOther.asSamTag());

        read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, saJoiner.toString());
        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));

        // Note that when we are unmapping a supplementary we do not bother unsetting all of its properties, because we will drop an
        // unmapped supplementary read immediately.
        assertTrue(checkTransformRead(read, CHR_3));
        assertTrue(read.getReadUnmappedFlag());
    }

    @Test
    public void testUnmapSupplementaryWhenPrimaryIsUnmappedUnpaired()
    {
        // Note that the suppReadData for a supplementary read refers to the associated primary read.
        SupplementaryReadData suppReadData = new SupplementaryReadData(CHR_3, 550, SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

        SAMRecord read = SamRecordTestUtils.createSamRecordUnpaired(
                TEST_READ_ID, CHR_2, 100, READ_BASES, TEST_READ_CIGAR, false, true, suppReadData);

        // Note that when we are unmapping a supplementary we do not bother unsetting all of its properties, because we will drop an unmapped
        // supplementary read immediately.
        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(checkTransformRead(read, CHR_2));
        assertTrue(read.getReadUnmappedFlag());
    }
}
