package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createRegion;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findUniqueBases;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class ReadCountsTest
{
    public static final String REF_BASE_STR_1 = "ABCDEFGHIJKLMNOPQRST";
    public static final String REF_BASE_STR_2 = REF_BASE_STR_1 + REF_BASE_STR_1;

    @Test
    public void testReadCoordinates()
    {
        // single matched region
        ReadRecord read = createReadRecord(1, "1", 100, 119, REF_BASE_STR_1, createCigar(0, 20, 0));

        // simple matched read
        assertEquals(1, read.getMappedRegionCoords().size());
        assertEquals(100, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(119, read.getMappedRegionCoords().get(0)[SE_END]);

        // with soft-clippings
        read = createReadRecord(1, "1", 100, 119, REF_BASE_STR_1, createCigar(10, 20, 10));

        assertEquals(1, read.getMappedRegionCoords().size());
        assertEquals(100, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(119, read.getMappedRegionCoords().get(0)[SE_END]);

        // with 2 splits
        Cigar cigar = new Cigar();
        cigar.add(new CigarElement(5, CigarOperator.M)); // matches half of first exon
        cigar.add(new CigarElement(20, CigarOperator.N));
        cigar.add(new CigarElement(10, CigarOperator.M)); // matches all of second exon
        cigar.add(new CigarElement(30, CigarOperator.N));
        cigar.add(new CigarElement(5, CigarOperator.M)); // matches past last exon into intron

        read = createReadRecord(1, "1", 100, 159, REF_BASE_STR_1, cigar);

        assertEquals(3, read.getMappedRegionCoords().size());
        assertEquals(100, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(104, read.getMappedRegionCoords().get(0)[SE_END]);
        assertEquals(125, read.getMappedRegionCoords().get(1)[SE_START]);
        assertEquals(134, read.getMappedRegionCoords().get(1)[SE_END]);
        assertEquals(165, read.getMappedRegionCoords().get(2)[SE_START]);
        assertEquals(169, read.getMappedRegionCoords().get(2)[SE_END]);

        // with a delete
        // 10M2D8M

        cigar = new Cigar();
        cigar.add(new CigarElement(10, CigarOperator.M));
        cigar.add(new CigarElement(2, CigarOperator.D));
        cigar.add(new CigarElement(8, CigarOperator.M));

        read = createReadRecord(1, "1", 100, 119, REF_BASE_STR_1.substring(0, 18), cigar);

        assertEquals(1, read.getMappedRegionCoords().size());
        assertEquals(100, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(119, read.getMappedRegionCoords().get(0)[SE_END]);

        // with an insert
        // 10M2D8M

        cigar = new Cigar();
        cigar.add(new CigarElement(10, CigarOperator.M));
        cigar.add(new CigarElement(2, CigarOperator.I));
        cigar.add(new CigarElement(8, CigarOperator.M));

        read = createReadRecord(1, "1", 100, 117, REF_BASE_STR_1, cigar);

        assertEquals(1, read.getMappedRegionCoords().size());
        assertEquals(100, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(117, read.getMappedRegionCoords().get(0)[SE_END]);
    }

    @Test
    public void testReadRegionTypes()
    {
        String REF_BASE_STR_1 = "ABCDEFGHIJKLMNOPQRST"; // currently unused

        // single read and exon
        RegionReadData region = new RegionReadData("1", 100, 200);

        // read covers part of the exon
        ReadRecord read = createReadRecord(1, "1", 110, 130, REF_BASE_STR_1, createCigar(0, 21, 0));

        assertEquals(1, read.getMappedRegionCoords().size());
        assertEquals(110, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(130, read.getMappedRegionCoords().get(0)[SE_END]);

        // test classification of reads
        assertEquals(WITHIN_EXON, read.getRegionMatchType(region));

        read = createReadRecord(1, "1", 90, 110, REF_BASE_STR_1, createCigar(0, 21, 0));
        assertEquals(EXON_INTRON, read.getRegionMatchType(region));

        read = createReadRecord(1, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 101, 0));
        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region));

        read = createReadRecord(1, "1", 100, 150, REF_BASE_STR_1, createCigar(0, 51, 0));
        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region));

        // read covering multiple exons
        Cigar cigar = new Cigar();
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches half of first exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(21, CigarOperator.M)); // matches all of second exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(25, CigarOperator.M)); // matches past last exon into intron

        read = createReadRecord(1, "1", 110, 204, REF_BASE_STR_1, cigar);

        assertEquals(3, read.getMappedRegionCoords().size());
        assertEquals(110, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(120, read.getMappedRegionCoords().get(0)[SE_END]);
        assertEquals(140, read.getMappedRegionCoords().get(1)[SE_START]);
        assertEquals(160, read.getMappedRegionCoords().get(1)[SE_END]);
        assertEquals(180, read.getMappedRegionCoords().get(2)[SE_START]);
        assertEquals(204, read.getMappedRegionCoords().get(2)[SE_END]);

        RegionReadData region1 = new RegionReadData("1", 100, 120);
        RegionReadData region2 = new RegionReadData("1", 140, 160);
        RegionReadData region3 = new RegionReadData("1", 180, 200);

        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region1));
        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region2));
        assertEquals(EXON_INTRON, read.getRegionMatchType(region3));
    }

    @Test
    public void testUnmappedSpliceJunctions()
    {
        // soft-clipping at one end
        RegionReadData region1 = new RegionReadData("1", 141, 150);
        region1.setRefBases(REF_BASE_STR_1.substring(0, 10));
        RegionReadData region2 = new RegionReadData("1", 200, 209);
        region2.setRefBases(REF_BASE_STR_1.substring(10, 19));

        region1.addPostRegion(region2);
        region2.addPreRegion(region1);

        ReadRecord read = createReadRecord(1, "1", 200, 209, REF_BASE_STR_1.substring(5, 19),
                createCigar(5, 10, 0));

        List<RegionReadData> regions = Lists.newArrayList(region2);
        read.processOverlappingRegions(regions);

        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region1));
        assertEquals(2, read.getMappedRegionCoords().size());
        assertEquals(146, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(150, read.getMappedRegionCoords().get(0)[SE_END]);

        // test again on the up side, with 2 different regions matching the inferred bases
        read = createReadRecord(1, "1", 141, 155, REF_BASE_STR_1.substring(0, 15),
                createCigar(0, 15, 0));

        RegionReadData region3 = new RegionReadData("1", 200, 229);
        region3.setRefBases(REF_BASE_STR_1.substring(10, 19) + REF_BASE_STR_1);

        RegionReadData region4 = new RegionReadData("1", 121, 150);
        region4.setRefBases(REF_BASE_STR_1 + REF_BASE_STR_1.substring(0, 10));

        region1.addPostRegion(region3);
        region3.addPreRegion(region1);

        region4.addPostRegion(region2);
        region4.addPostRegion(region3);
        region2.addPreRegion(region4);
        region3.addPreRegion(region4);

        regions = Lists.newArrayList(region1, region4);
        read.processOverlappingRegions(regions);

        assertEquals(EXON_BOUNDARY, read.getMappedRegions().get(region1));
        assertEquals(EXON_BOUNDARY, read.getMappedRegions().get(region2));
        assertEquals(EXON_BOUNDARY, read.getMappedRegions().get(region3));
        assertEquals(EXON_BOUNDARY, read.getMappedRegions().get(region4));
        assertEquals(2, read.getMappedRegionCoords().size());
        assertEquals(200, read.getMappedRegionCoords().get(read.getMappedRegionCoords().size() - 1)[SE_START]);
        assertEquals(204, read.getMappedRegionCoords().get(read.getMappedRegionCoords().size() - 1)[SE_END]);

        // test with a soft-clipped and overhanging base together
        read = createReadRecord(1, "1", 141, 152, REF_BASE_STR_1.substring(0, 15),
                createCigar(0, 12, 3));

        read.processOverlappingRegions(regions);

        assertEquals(EXON_BOUNDARY, read.getMappedRegions().get(region1));
        assertEquals(EXON_BOUNDARY, read.getMappedRegions().get(region2));
        assertEquals(EXON_BOUNDARY, read.getMappedRegions().get(region3));
        assertEquals(EXON_BOUNDARY, read.getMappedRegions().get(region4));
        assertEquals(2, read.getMappedRegionCoords().size());
        assertEquals(200, read.getMappedRegionCoords().get(read.getMappedRegionCoords().size() - 1)[SE_START]);
        assertEquals(204, read.getMappedRegionCoords().get(read.getMappedRegionCoords().size() - 1)[SE_END]);

        // test again with ambiguous mapping of 1 base to more than 1 adjacent exons, and observe the read coords being truncated
        read = createReadRecord(1, "1", 132, 151,
                REF_BASE_STR_1.substring(0, 19) + REF_BASE_STR_1.substring(10, 11), createCigar(0, 20, 0));

        read.processOverlappingRegions(Lists.newArrayList(region1, region4));
        assertEquals(1, read.getMappedRegionCoords().size());
        assertEquals(150, read.getMappedRegionCoords().get(0)[SE_END]);

        read = createReadRecord(1, "1", 199, 218,
                REF_BASE_STR_1.substring(9, 10) + REF_BASE_STR_1.substring(0, 19), createCigar(0, 20, 0));

        read.processOverlappingRegions(Lists.newArrayList(region2, region3));
        assertEquals(1, read.getMappedRegionCoords().size());
        assertEquals(200, read.getMappedRegionCoords().get(0)[SE_START]);
    }

    @Test
    public void testUniqueRegionBases()
    {
        RegionReadData region1 = createRegion("GEN01", 1, 1, "1", 100, 119);
        region1.setRefBases(REF_BASE_STR_1);

        // covers the first region entirely
        RegionReadData region2 = createRegion("GEN01", 2, 1, "1", 95, 124);
        region2.setRefBases(REF_BASE_STR_1 + REF_BASE_STR_1.substring(0, 10));

        RegionReadData region3 = createRegion("GEN01", 3, 1, "1", 80, 99);
        region3.setRefBases(REF_BASE_STR_1);

        RegionReadData region4 = createRegion("GEN01", 4, 1, "1", 70, 89);
        region4.setRefBases(REF_BASE_STR_1);

        RegionReadData region5 = createRegion("GEN01", 5, 1, "1", 130, 149);
        region5.setRefBases(REF_BASE_STR_1);

        List<RegionReadData> regions = Lists.newArrayList(region1, region2, region3, region4, region5);
        findUniqueBases(regions);

        assertEquals(0, region1.uniqueBaseCount());
        assertEquals(5, region2.uniqueBaseCount());
        assertEquals(5, region3.uniqueBaseCount());
        assertEquals(10, region4.uniqueBaseCount());
        assertEquals(20, region5.uniqueBaseCount());
    }

}