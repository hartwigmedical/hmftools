package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.markdups.common.ReadUnmapper;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class UnmapReadsTest
{
    @Test
    public void testUnamppingCoords()
    {
        Map<String,List<BaseRegion>> chrLocationsMap = Maps.newHashMap();
        chrLocationsMap.put(CHR_1, Lists.newArrayList(new BaseRegion(500, 700)));

        ReadUnmapper readUnmapper = new ReadUnmapper(chrLocationsMap, 100);

        String readBases = "";
        String readCigar = "100M";
        String readId = "READ_001";

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 100, readBases, readCigar, CHR_2, 600, false,
                false, null, false, readCigar);

        assertFalse(readUnmapper.checkTransformRead(read, chrLocationsMap.get(CHR_1)));

        // read in region
        read = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 550, readBases, readCigar, CHR_2, 600, false,
                false, null, false, readCigar);

        assertTrue(readUnmapper.checkTransformRead(read, chrLocationsMap.get(CHR_1)));

        assertTrue(read.getReadUnmappedFlag());
        assertEquals(600, read.getAlignmentStart());
        assertEquals(CHR_2, read.getReferenceName());

        // mate in region
        read = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 100, readBases, readCigar, CHR_1, 600, false,
                false, null, false, readCigar);

        assertTrue(readUnmapper.checkTransformRead(read, chrLocationsMap.get(CHR_1)));

        assertFalse(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
        assertEquals(100, read.getMateAlignmentStart());
        assertEquals(CHR_1, read.getMateReferenceName());

        // drop supp alignment
        read = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 100, readBases, readCigar, CHR_1, 600, false,
                false,
                new SupplementaryReadData(CHR_1, 550, SUPP_POS_STRAND, readCigar, 60),
                false, readCigar);

        assertTrue(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));

        assertTrue(readUnmapper.checkTransformRead(read, chrLocationsMap.get(CHR_1)));
        assertFalse(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE));

        // read goes to fully unmapped
        read = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 550, readBases, readCigar, NO_CHROMOSOME_NAME, NO_POSITION, false,
                false, null, false, readCigar);

        assertTrue(readUnmapper.checkTransformRead(read, chrLocationsMap.get(CHR_1)));

        assertTrue(read.getReadUnmappedFlag());
        assertTrue(read.getMateUnmappedFlag());
    }

    @Test
    public void testUnamppingScenarios()
    {


    }

}
