package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createRegion;
import static com.hartwig.hmftools.isofox.common.ReadRecord.markRegionBases;
import static com.hartwig.hmftools.isofox.common.CommonUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.CommonUtils.findStringOverlaps;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;

import org.junit.Test;

import htsjdk.samtools.Cigar;

public class ReadUtilsTest
{
    @Test
    public void testCigarCreation()
    {
        Cigar cigar = createCigar(2, 10, 1);
        assertTrue(cigar.toString().equals("2S10M1S"));

        cigar = createCigar(0, 10, 100, 12, 0);
        assertTrue(cigar.toString().equals("10M100N12M"));

        cigar = createCigar(2, 10, 100, 12, 4);
        assertTrue(cigar.toString().equals("2S10M100N12M4S"));
    }

    @Test
    public void testMappingCoords()
    {
        List<int[]> mappings1 = Lists.newArrayList();

        // no overlaps
        mappings1.add(new int[] { 10, 20 });
        mappings1.add(new int[] { 40, 50 });

        List<int[]> mappings2 = Lists.newArrayList();

        mappings2.add(new int[] { 60, 70 });
        mappings2.add(new int[] { 80, 90 });

        List<int[]> commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(4, commonMappings.size());

        mappings1.clear();
        mappings2.clear();

        // widening of all regions only
        mappings1.add(new int[] { 10, 20 });
        mappings1.add(new int[] { 40, 50 });
        mappings1.add(new int[] { 70, 80 });

        // no overlaps
        mappings2.add(new int[] { 25, 35 });
        mappings2.add(new int[] { 55, 65 });
        mappings2.add(new int[] { 85, 95 });

        commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(6, commonMappings.size());
        assertTrue(commonMappings.contains(mappings1.get(0)));
        assertTrue(commonMappings.contains(mappings1.get(1)));
        assertTrue(commonMappings.contains(mappings1.get(2)));
        assertTrue(commonMappings.contains(mappings2.get(0)));
        assertTrue(commonMappings.contains(mappings2.get(1)));
        assertTrue(commonMappings.contains(mappings2.get(2)));

        // widening of all regions only
        mappings2.clear();

        mappings2.add(new int[] { 5, 15 });
        mappings2.add(new int[] { 35, 45 });
        mappings2.add(new int[] { 55, 75 });

        commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(3, commonMappings.size());
        assertEquals(5, commonMappings.get(0)[SE_START]);
        assertEquals(20, commonMappings.get(0)[SE_END]);
        assertEquals(35, commonMappings.get(1)[SE_START]);
        assertEquals(50, commonMappings.get(1)[SE_END]);
        assertEquals(55, commonMappings.get(2)[SE_START]);
        assertEquals(80, commonMappings.get(2)[SE_END]);

        // one other region overlapping all others
        mappings2.clear();

        mappings2.add(new int[] { 5, 95 });

        commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(1, commonMappings.size());
        assertEquals(5, commonMappings.get(0)[SE_START]);
        assertEquals(95, commonMappings.get(0)[SE_END]);

        mappings2.clear();
        mappings1.clear();

        // a mix of various scenarios
        mappings1.add(new int[] { 10, 20 });

        mappings2.add(new int[] { 30, 40 });

        mappings2.add(new int[] { 50, 60 });
        mappings1.add(new int[] { 55, 75 });
        mappings1.add(new int[] { 85, 95 });
        mappings2.add(new int[] { 70, 110 });

        mappings2.add(new int[] { 120, 130 });

        mappings1.add(new int[] { 140, 150 });

        commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(5, commonMappings.size());

        assertEquals(50, commonMappings.get(2)[SE_START]);
        assertEquals(110, commonMappings.get(2)[SE_END]);
    }

    @Test
    public void testFragmentTracking()
    {
        FragmentTracker fragTracker = new FragmentTracker();

        String readId1 = "read1";
        String readId2 = "read2";
        String readId3 = "read3";

        assertFalse(fragTracker.checkReadId(readId1));
        assertFalse(fragTracker.checkReadId(readId2));
        assertFalse(fragTracker.checkReadId(readId3));

        assertEquals(3, fragTracker.readsCount());

        assertTrue(fragTracker.checkReadId(readId1));
        assertEquals(2, fragTracker.readsCount());
        assertTrue(fragTracker.checkReadId(readId2));
        assertEquals(1, fragTracker.readsCount());
        assertTrue(fragTracker.checkReadId(readId3));
        assertEquals(0, fragTracker.readsCount());

        ReadRecord read1 = createReadRecord(1, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));
        ReadRecord read2 = createReadRecord(2, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));
        ReadRecord read3 = createReadRecord(3, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));

        assertEquals(null, fragTracker.checkRead(read1));
        assertEquals(null, fragTracker.checkRead(read2));
        assertEquals(null, fragTracker.checkRead(read3));

        assertEquals(3, fragTracker.readsCount());

        ReadRecord read1b = createReadRecord(1, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));
        ReadRecord read2b = createReadRecord(2, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));
        ReadRecord read3b = createReadRecord(3, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));

        assertEquals(read1, fragTracker.checkRead(read1b));
        assertEquals(read2, fragTracker.checkRead(read2b));
        assertEquals(read3, fragTracker.checkRead(read3b));

        assertEquals(0, fragTracker.readsCount());
    }

    @Test
    public void testBaseAssignment()
    {
        RegionReadData region = createRegion("GEN01", 1, 1, "1", 100, 119);
        region.setRefBases(REF_BASE_STR_1);

        List<int[]> readCoords = Lists.newArrayList();
        readCoords.add(new int[] { 100, 119 });

        markRegionBases(readCoords, region);
        assertEquals(20, region.baseCoverage(1));

        region.clearState();

        readCoords.clear();
        readCoords.add(new int[] { 100, 104 });
        readCoords.add(new int[] { 110, 114 });
        readCoords.add(new int[] { 118, 119 });

        markRegionBases(readCoords, region);
        assertEquals(12, region.baseCoverage(1));
    }

    @Test
    public void testBaseComparisons()
    {
        // extra bases at the start
        String str1 = "ABCDEFGHIJ";
        String str2 = "XXXABCDEFGHIJ";

        int overlap = findStringOverlaps(str1, str2);
        assertEquals(10, overlap);

        overlap = findStringOverlaps(str2, str1);
        assertEquals(10, overlap);

        // and in the middle
        str1 = "ABCDEFZZZGHIJ";
        str2 = "XXXABCDEFGHIJ";

        overlap = findStringOverlaps(str1, str2);
        assertEquals(10, overlap);

        overlap = findStringOverlaps(str2, str1);
        assertEquals(10, overlap);

        // some incorrect letters - 2/21 is more than 90%
        str1 = "ABCDEFGHIYKLMNOPQRSTU";
        str2 = "ABCXEFGHIJKLMNOPQRSTU";

        overlap = findStringOverlaps(str1, str2);
        assertEquals(19, overlap);
    }
}
