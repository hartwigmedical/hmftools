package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_2;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createRegion;
import static com.hartwig.hmftools.isofox.common.Read.markRegionBases;
import static com.hartwig.hmftools.isofox.common.CommonUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.CommonUtils.findStringOverlaps;
import static com.hartwig.hmftools.isofox.common.ReadUtils.trimAdapterBases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.MappedCoords;
import com.hartwig.hmftools.isofox.common.Read;
import com.hartwig.hmftools.isofox.common.RegionReadData;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFlag;

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
        List<CigarElement> cigarElements = Lists.newArrayList(
                new CigarElement(10, M),
                new CigarElement(10, I),
                new CigarElement(10, M),
                new CigarElement(10, D),
                new CigarElement(10, M));

        MappedCoords mappedCoords = MappedCoords.build(cigarElements, 100);
        assertEquals(100, mappedCoords.alignments().get(0).start());
        assertEquals(139, mappedCoords.alignments().get(0).end());

        assertEquals(1, mappedCoords.alignmentCount());

        assertTrue(mappedCoords.alignmentsOverlap(90, 101));
        assertTrue(mappedCoords.alignmentsOverlap(138, 150));
        assertTrue(mappedCoords.alignmentsWithin(99, 139));

        cigarElements = Lists.newArrayList(
                new CigarElement(11, M),
                new CigarElement(9, N),
                new CigarElement(11, M),
                new CigarElement(9, N),
                new CigarElement(11, M));

        mappedCoords = MappedCoords.build(cigarElements, 100);

        assertEquals(3, mappedCoords.alignmentCount());
        assertEquals(100, mappedCoords.alignments().get(0).start());
        assertEquals(110, mappedCoords.alignments().get(0).end());
        assertEquals(120, mappedCoords.alignments().get(1).start());
        assertEquals(130, mappedCoords.alignments().get(1).end());
        assertEquals(140, mappedCoords.alignments().get(2).start());
        assertEquals(150, mappedCoords.alignments().get(2).end());

        assertEquals(3, mappedCoords.alignmentCount());
        assertEquals(3, mappedCoords.originalAlignmentCount());

        // now add inferred sections
        mappedCoords.addInferredRegion(true, 80, 90);

        assertEquals(4, mappedCoords.alignmentCount());
        assertEquals(3, mappedCoords.originalAlignmentCount());

        mappedCoords.addInferredRegion(false, 160, 170);
        assertEquals(5, mappedCoords.alignmentCount());
    }

    @Test
    public void testOverlappingCoordinates()
    {
        List<BaseRegion> mappings1 = Lists.newArrayList();

        // no overlaps
        mappings1.add(new BaseRegion(10, 20));
        mappings1.add(new BaseRegion(40, 50));

        List<BaseRegion> mappings2 = Lists.newArrayList();

        mappings2.add(new BaseRegion(60, 70));
        mappings2.add(new BaseRegion(80, 90));

        List<BaseRegion> commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(4, commonMappings.size());

        mappings1.clear();
        mappings2.clear();

        // widening of all regions only
        mappings1.add(new BaseRegion(10, 20));
        mappings1.add(new BaseRegion(40, 50));
        mappings1.add(new BaseRegion(70, 80));

        // no overlaps
        mappings2.add(new BaseRegion(25, 35));
        mappings2.add(new BaseRegion(55, 65));
        mappings2.add(new BaseRegion(85, 95));

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

        mappings2.add(new BaseRegion(5, 15));
        mappings2.add(new BaseRegion(35, 45));
        mappings2.add(new BaseRegion(55, 75));

        commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(3, commonMappings.size());
        assertEquals(5, commonMappings.get(0).start());
        assertEquals(20, commonMappings.get(0).end());
        assertEquals(35, commonMappings.get(1).start());
        assertEquals(50, commonMappings.get(1).end());
        assertEquals(55, commonMappings.get(2).start());
        assertEquals(80, commonMappings.get(2).end());

        // one other region overlapping all others
        mappings2.clear();

        mappings2.add(new BaseRegion(5, 95));

        commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(1, commonMappings.size());
        assertEquals(5, commonMappings.get(0).start());
        assertEquals(95, commonMappings.get(0).end());

        mappings2.clear();
        mappings1.clear();

        // a mix of various scenarios
        mappings1.add(new BaseRegion(10, 20));

        mappings2.add(new BaseRegion(30, 40));

        mappings2.add(new BaseRegion(50, 60));
        mappings1.add(new BaseRegion(55, 75));
        mappings1.add(new BaseRegion(85, 95));
        mappings2.add(new BaseRegion(70, 110));

        mappings2.add(new BaseRegion(120, 130));

        mappings1.add(new BaseRegion(140, 150));

        commonMappings = deriveCommonRegions(mappings1, mappings2);
        assertEquals(5, commonMappings.size());

        assertEquals(50, commonMappings.get(2).start());
        assertEquals(110, commonMappings.get(2).end());
    }

    @Test
    public void testFragmentTracking()
    {
        FragmentTracker fragTracker = new FragmentTracker();

        Read read1 = createReadRecord(1, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));
        Read read2 = createReadRecord(2, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));
        Read read3 = createReadRecord(3, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));

        assertEquals(null, fragTracker.checkRead(read1));
        assertEquals(null, fragTracker.checkRead(read2));
        assertEquals(null, fragTracker.checkRead(read3));

        assertEquals(3, fragTracker.readsCount());

        Read read1b = createReadRecord(1, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));
        Read read2b = createReadRecord(2, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));
        Read read3b = createReadRecord(3, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 50, 0));

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

        List<BaseRegion> readCoords = Lists.newArrayList();
        readCoords.add(new BaseRegion(100, 119));

        markRegionBases(readCoords, region);
        assertEquals(20, region.baseCoverage(1));

        region.clearState();

        readCoords.clear();
        readCoords.add(new BaseRegion(100, 104));
        readCoords.add(new BaseRegion(110, 114));
        readCoords.add(new BaseRegion(118, 119));

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

    @Test
    public void testAdapterTrimming2()
    {
        Read read1 = createReadRecord(1, CHR_1, 105, 124, REF_BASE_STR_2, createCigar(10, 20, 10));
        Read read2 = createReadRecord(1, CHR_1, 100, 119, REF_BASE_STR_2, createCigar(10, 20, 10));
        read2.setFlag(SAMFlag.READ_REVERSE_STRAND, true);

        assertEquals(95, read1.unclippedStart());
        assertEquals(134, read1.unclippedEnd());
        assertEquals(90, read2.unclippedStart());
        assertEquals(129, read2.unclippedEnd());

        trimAdapterBases(read1, read2);

        assertEquals(129, read1.unclippedEnd());
        assertEquals(95, read2.unclippedStart());

        assertEquals("10S20M5S", read1.cigarStr());
        assertEquals("5S20M10S", read2.cigarStr());

        // test passing in reversed
        read1 = createReadRecord(1, CHR_1, 105, 124, REF_BASE_STR_2, createCigar(10, 20, 10));
        read2 = createReadRecord(1, CHR_1, 100, 119, REF_BASE_STR_2, createCigar(10, 20, 10));
        read2.setFlag(SAMFlag.READ_REVERSE_STRAND, true);

        trimAdapterBases(read2, read1);

        assertEquals(129, read1.unclippedEnd());
        assertEquals(95, read2.unclippedStart());

        assertEquals("10S20M5S", read1.cigarStr());
        assertEquals("5S20M10S", read2.cigarStr());

        // cap trimming at soft-clips
        read1 = createReadRecord(1, CHR_1, 105, 124, REF_BASE_STR_2, createCigar(10, 20, 10));
        read2 = createReadRecord(1, CHR_1, 94, 113, REF_BASE_STR_2, createCigar(10, 20, 10));
        read2.setFlag(SAMFlag.READ_REVERSE_STRAND, true);

        assertEquals(95, read1.unclippedStart());
        assertEquals(134, read1.unclippedEnd());
        assertEquals(84, read2.unclippedStart());
        assertEquals(123, read2.unclippedEnd());

        trimAdapterBases(read1, read2);

        assertEquals(124, read1.unclippedEnd());
        assertEquals(94, read2.unclippedStart());

        assertEquals("10S20M", read1.cigarStr());
        assertEquals("20M10S", read2.cigarStr());

        // test a scenario where N-splitting has lead to conflicting alignments, requiring position indices to find the point of adapter SCs

        // test 1: no trimming
        // read 1: 100-114 - 886N - 1000-1019 - 5S
        // read 2: 5S -             1000-1019 - 981N - 2000-2014

        read1 = createReadRecord(1, CHR_1, 100, 1019, REF_BASE_STR_2, cigarFromStr("15M885N20M5S"));
        read2 = createReadRecord(1, CHR_1, 1000, 2014, REF_BASE_STR_2, cigarFromStr("5S20M980N15M"));
        read2.setFlag(SAMFlag.READ_REVERSE_STRAND, true);

        assertEquals(100, read1.unclippedStart());
        assertEquals(1024, read1.unclippedEnd());
        assertEquals(995, read2.unclippedStart());
        assertEquals(2014, read2.unclippedEnd());

        trimAdapterBases(read2, read1);

        assertEquals(100, read1.unclippedStart());
        assertEquals(1024, read1.unclippedEnd());
        assertEquals(995, read2.unclippedStart());
        assertEquals(2014, read2.unclippedEnd());

        // test 2: trim from both sides
        // read 1: 110-114 - 886N - 1000-1019 - 10S
        // read 2: 10S -             1000-1019 - 981N - 2000-2004
        read1 = createReadRecord(1, CHR_1, 110, 1019, REF_BASE_STR_2, cigarFromStr("5M885N20M10S"));
        read2 = createReadRecord(1, CHR_1, 1000, 2004, REF_BASE_STR_2, cigarFromStr("10S20M980N5M"));
        read2.setFlag(SAMFlag.READ_REVERSE_STRAND, true);

        assertEquals(110, read1.unclippedStart());
        assertEquals(1029, read1.unclippedEnd());
        assertEquals(990, read2.unclippedStart());
        assertEquals(2004, read2.unclippedEnd());

        trimAdapterBases(read2, read1);

        assertEquals(1024, read1.unclippedEnd());
        assertEquals(995, read2.unclippedStart());

        assertEquals("5M885N20M5S", read1.cigarStr());
        assertEquals("5S20M980N5M", read2.cigarStr());
    }
}