package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.markdups.common.Constants.UNMAP_MIN_HIGH_DEPTH;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.markdups.common.HighDepthRegion;
import com.hartwig.hmftools.markdups.common.ReadUnmapper;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UnmapReadsTest
{
    private static final String READ_BASES = "A".repeat(100);
    private static final String READ_CIGAR = "100M";
    private static final String SOFT_CLIPPED_READ_CIGAR = "21S79M";
    private static final String READ_ID = "READ_001";
    private static final Map<String, List<HighDepthRegion>> CHR_LOCATION_MAP;
    private static final ReadUnmapper READ_UNMAPPER;

    static
    {
        CHR_LOCATION_MAP = Maps.newHashMap();
        CHR_LOCATION_MAP.put(CHR_1, Lists.newArrayList(new HighDepthRegion(500, 700, 0)));
        CHR_LOCATION_MAP.put(CHR_3, Lists.newArrayList(new HighDepthRegion(500, 700, UNMAP_MIN_HIGH_DEPTH)));
        READ_UNMAPPER = new ReadUnmapper(CHR_LOCATION_MAP);
    }

    @Test
    public void testGetSoftClipCountFromCigarStr()
    {
        assertEquals(0, ReadUnmapper.getSoftClipCountFromCigarStr("151M"));
        assertEquals(10, ReadUnmapper.getSoftClipCountFromCigarStr("10S141M"));
        assertEquals(10, ReadUnmapper.getSoftClipCountFromCigarStr("141M10S"));
        assertEquals(12, ReadUnmapper.getSoftClipCountFromCigarStr("2S138M10S"));
    }

    @Test
    public void testNoUnmap()
    {
        // neither read in unmapp regions
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_2, 600, false,
                false, null, true, READ_CIGAR);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertFalse(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertFalse(read.getReadUnmappedFlag());
        assertFalse(read.getMateUnmappedFlag());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testNoUnmapSoftClipped()
    {
        // soft-clipped but not in an unmap region
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 100, READ_BASES, SOFT_CLIPPED_READ_CIGAR, CHR_2, 600, false,
                false, null, true, SOFT_CLIPPED_READ_CIGAR);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertFalse(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertFalse(read.getReadUnmappedFlag());
        assertFalse(read.getMateUnmappedFlag());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapReadSoftClipped()
    {
        // read unmapped since chimeric, sets read coords to mate
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 550, READ_BASES, SOFT_CLIPPED_READ_CIGAR, CHR_2, 600, false,
                false, null, true, READ_CIGAR);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertTrue(read.getReadUnmappedFlag());
        assertFalse(read.getMateUnmappedFlag());
        assertEquals(600, read.getAlignmentStart());
        assertEquals(CHR_2, read.getReferenceName());
    }

    @Test
    public void testUnmapMateSoftClipped()
    {
        // read unmapped since soft-cliped
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_1, 600, false,
                false, null, true, SOFT_CLIPPED_READ_CIGAR);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertFalse(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(100, read.getMateAlignmentStart());
        assertEquals(CHR_1, read.getMateReferenceName());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapMateAndSupplementarySoftClipped()
    {
        // mate unmapped and supp removed since both are in unmap regions and both have long SCs
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_1, 600, false,
                false,
                new SupplementaryReadData(CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60),
                true, SOFT_CLIPPED_READ_CIGAR);

        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertFalse(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(100, read.getMateAlignmentStart());
        assertEquals(CHR_1, read.getMateReferenceName());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapReadSoftClippedMateUnmapped()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 550, READ_BASES, SOFT_CLIPPED_READ_CIGAR, CHR_1, 550, false,
                false, null, true, null);
        read.setMateUnmappedFlag(true);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertTrue(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(NO_POSITION, read.getAlignmentStart());
        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
        assertEquals(NO_POSITION, read.getMateAlignmentStart());
        assertEquals(NO_CHROMOSOME_NAME, read.getMateReferenceName());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapMateSoftClippedReadUnmapped()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 550, READ_BASES, NO_CIGAR, CHR_1, 550, false,
                false, null, true, SOFT_CLIPPED_READ_CIGAR);
        read.setReadUnmappedFlag(true);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertTrue(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(NO_POSITION, read.getAlignmentStart());
        assertEquals(NO_CHROMOSOME_NAME, read.getReferenceName());
        assertEquals(NO_POSITION, read.getMateAlignmentStart());
        assertEquals(NO_CHROMOSOME_NAME, read.getMateReferenceName());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapSupplementaryWhenPrimaryIsUnmappedSoftClipped()
    {
        // Note that the suppReadData for a supplementary read refers to the associated primary read
        SupplementaryReadData suppReadData = new SupplementaryReadData(
                CHR_1, 550, SUPP_POS_STRAND, SOFT_CLIPPED_READ_CIGAR, 60);

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_2, 100, READ_BASES, READ_CIGAR, CHR_2, 1000, false,
                true, suppReadData, true, READ_CIGAR);

        // when supplementaries are unmapped not all properties are unset since it will be dropped immediately
        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertTrue(read.getReadUnmappedFlag());
    }

    @Test
    public void testUnmapReadChimeric()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 550, READ_BASES, READ_CIGAR, CHR_2, 600, false,
                false, null, true, READ_CIGAR);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertTrue(read.getReadUnmappedFlag());
        assertFalse(read.getMateUnmappedFlag());
        assertEquals(600, read.getAlignmentStart());
        assertEquals(CHR_2, read.getReferenceName());
    }

    @Test
    public void testUnmapMateChimeric()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_2, 100, READ_BASES, READ_CIGAR, CHR_1, 600, false,
                false, null, true, READ_BASES);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertFalse(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(100, read.getMateAlignmentStart());
        assertEquals(CHR_2, read.getMateReferenceName());
    }

    @Test
    public void testSupplementaryNotUnmappedChimeric()
    {
        // Note that the suppReadData for a supplementary read refers to the associated primary read.
        SupplementaryReadData suppReadData = new SupplementaryReadData(CHR_2, 100, SUPP_POS_STRAND, READ_CIGAR, 60);

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_1, 100, READ_BASES, READ_CIGAR, CHR_2, 600, false, true,
                suppReadData, true, READ_BASES);

        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertFalse(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertFalse(read.getReadUnmappedFlag());
        assertFalse(read.getMateUnmappedFlag());
        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }

    @Test
    public void testUnmapSupplementaryWhenPrimaryIsUnmappedChimeric()
    {
        // Note that the suppReadData for a supplementary read refers to the associated primary read.
        SupplementaryReadData suppReadData = new SupplementaryReadData(CHR_1, 550, SUPP_POS_STRAND, READ_CIGAR, 60);

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_2, 100, READ_BASES, READ_CIGAR, CHR_2, 1000, false,
                true, suppReadData, true, READ_CIGAR);

        // Note that when we are unmapping a supplementary we do not bother unsetting all of its properties, because we will drop an unmapped
        // supplementary read immediately.
        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_1)));
        assertTrue(read.getReadUnmappedFlag());
    }

    @Test
    public void testUnmapReadHighDepth()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_3, 600, READ_BASES, READ_CIGAR, CHR_3, 701, false,
                false, null, true, READ_CIGAR);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_3)));
        assertTrue(read.getReadUnmappedFlag());
        assertFalse(read.getMateUnmappedFlag());
        assertEquals(701, read.getAlignmentStart());
        assertEquals(CHR_3, read.getReferenceName());
    }

    @Test
    public void testUnmapMateHighDepth()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_3, 400, READ_BASES, READ_CIGAR, CHR_3, 500, false,
                false, null, true, READ_CIGAR);

        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_3)));
        assertFalse(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(400, read.getMateAlignmentStart());
        assertEquals(CHR_3, read.getMateReferenceName());
    }

    @Test
    public void testUnmapMateAndSupplementaryHighDepth()
    {
        SAMRecord read = SamRecordTestUtils.createSamRecord(
                READ_ID, CHR_3, 400, READ_BASES, READ_CIGAR, CHR_3, 500, false,
                false,
                new SupplementaryReadData(CHR_3, 600, SUPP_POS_STRAND, READ_CIGAR, 60),
                true, READ_CIGAR);

        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
        assertTrue(READ_UNMAPPER.checkTransformRead(read, CHR_LOCATION_MAP.get(CHR_3)));
        assertFalse(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(400, read.getMateAlignmentStart());
        assertEquals(CHR_3, read.getMateReferenceName());
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));
    }
}
