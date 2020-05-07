package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneReadData;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createRegion;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.ReadRecord.markRegionBases;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_MATCH;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findUniqueBases;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.RnaUtils.findStringOverlaps;
import static com.hartwig.hmftools.isofox.common.TransMatchType.ALT;
import static com.hartwig.hmftools.isofox.common.TransMatchType.EXONIC;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.common.TransMatchType.UNSPLICED;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.SAMFlag.PROPER_PAIR;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.FragmentType;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;

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
        assertEquals(EXON_MATCH, read.getRegionMatchType(region));

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
        assertEquals(EXON_MATCH, read.getRegionMatchType(region2));
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

        ReadRecord read = createReadRecord(1, "1", 200, 209, REF_BASE_STR_1.substring(5, 19), createCigar(5, 10, 0));

        List<RegionReadData> regions = Lists.newArrayList(region2);
        read.processOverlappingRegions(regions);

        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region1));
        assertEquals(2, read.getMappedRegionCoords().size());
        assertEquals(146, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(150, read.getMappedRegionCoords().get(0)[SE_END]);

        // test again on the up side, with 2 different regions matching the inferred bases
        read = createReadRecord(1, "1", 141, 155, REF_BASE_STR_1.substring(0, 15), createCigar(0, 15, 0));

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
    public void testReadTranscriptClassification()
    {
        int trans1 = 1;

        RegionReadData region = createRegion("GEN01", trans1, 1, "1", 100, 200);

        // unspliced read does not support the transcript
        ReadRecord read = createReadRecord(1, "1", 90, 110, REF_BASE_STR_1, createCigar(0, 21, 0));

        List<RegionReadData> regions = Lists.newArrayList(region);
        read.processOverlappingRegions(regions);

        assertEquals(UNSPLICED, read.getTranscriptClassification(trans1));

        // exonic read does support the transcript
        read = createReadRecord(1, "1", 120, 140, REF_BASE_STR_1, createCigar(0, 21, 0));

        regions = Lists.newArrayList(region);
        read.processOverlappingRegions(regions);

        assertEquals(EXONIC, read.getTranscriptClassification(trans1));

        // skipped exon
        RegionReadData region1 = createRegion("GEN01", trans1, 1, "1", 100, 120);
        RegionReadData region2 = createRegion("GEN01", trans1, 2, "1", 140, 160);
        RegionReadData region3 = createRegion("GEN01", trans1, 3, "1", 180, 200);

        // read covering multiple exons but skips the middle exon
        read = createReadRecord(1, "1", 110, 200, REF_BASE_STR_1, createCigar(0, 11, 59, 21, 0));

        regions = Lists.newArrayList(region1, region2, region3);
        read.processOverlappingRegions(regions);

        assertEquals(ALT, read.getTranscriptClassification(trans1));

        // incomplete intermediary exon
        Cigar cigar = new Cigar();
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches half of first exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches part of second exon
        cigar.add(new CigarElement(29, CigarOperator.N));
        cigar.add(new CigarElement(25, CigarOperator.M)); // matches past last exon into intron

        read = createReadRecord(1, "1", 110, 200, REF_BASE_STR_1, cigar);

        read.processOverlappingRegions(regions);

        assertEquals(ALT, read.getTranscriptClassification(trans1));

        // read spanning into intron (unspliced)
        cigar = new Cigar();
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches half of first exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(31, CigarOperator.M)); // goes past second exon into intron
        cigar.add(new CigarElement(9, CigarOperator.N));
        cigar.add(new CigarElement(25, CigarOperator.M)); // matches past last exon into intron

        read = createReadRecord(1, "1", 110, 200, REF_BASE_STR_1, cigar);

        read.processOverlappingRegions(regions);

        assertEquals(UNSPLICED, read.getTranscriptClassification(trans1));

        // a read mapping to 2 regions, one as a splice junction, the other as exonic
        read = createReadRecord(1, "1", 181, 319, REF_BASE_STR_1 + REF_BASE_STR_1,
                createCigar(0, 20, 99, 20, 0));

        region1 = createRegion("GEN01", trans1, 1, "1", 100, 200);
        region2 = createRegion("GEN01", trans1, 2, "1", 300, 400);

        int trans2 = 2;
        region3 = createRegion("GEN01", trans2, 1, "1", 100, 220);

        regions = Lists.newArrayList(region1, region2, region3);
        read.processOverlappingRegions(regions);

        assertEquals(SPLICE_JUNCTION, read.getTranscriptClassification(trans1));
        assertEquals(ALT, read.getTranscriptClassification(trans2));
    }

    @Test
    public void testFragmentReadPairs()
    {
        // one read outside the gene
        IsofoxConfig config = new IsofoxConfig();
        BamFragmentAllocator bamReader = new BamFragmentAllocator(config, new ResultsWriter(config));

        int transId1 = 1;
        String transName1 = "TRANS01";
        String geneId = "GENE01";

        TranscriptData transData1 = new TranscriptData(transId1, transName1, geneId, true, (byte) 1,
                1000, 5000, null, null, "");

        transData1.exons().add(new ExonData(transId1, 1000, 1200, 1, -1, -1));
        transData1.exons().add(new ExonData(transId1, 2000, 2500, 2, -1, -1));
        transData1.exons().add(new ExonData(transId1, 4500, 5000, 3, -1, -1));

        GeneReadData geneReadData = createGeneReadData(geneId, "1", (byte) 1, 1000, 5000);
        geneReadData.setTranscripts(Lists.newArrayList(transData1));

        ReadRecord read1 = createReadRecord(1, "1", 100, 200, REF_BASE_STR_1, createCigar(0, 10, 0));
        ReadRecord read2 = createReadRecord(1, "1", 1050, 1150, REF_BASE_STR_1, createCigar(0, 20, 0));

        List<ReadRecord> reads = Lists.newArrayList(read1, read2);

        GeneCollection geneSet = new GeneCollection(0, Lists.newArrayList(geneReadData));
        bamReader.processReadRecords(geneSet, reads);

        int[] geneCounts = geneSet.getCounts();
        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(CHIMERIC)]);

        // exon to intronic read
        read1 = createReadRecord(1, "1", 1300, 1350, REF_BASE_STR_1, createCigar(0, 50, 0));
        read1.setFragmentInsertSize(1100);
        read2 = createReadRecord(1, "1", 2300, 2400, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-1100);

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(2, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(FragmentType.UNSPLICED)]);

        // fully intronic
        read1 = createReadRecord(1, "1", 2600, 2650, REF_BASE_STR_1, createCigar(0, 50, 0));
        read1.setFragmentInsertSize(300);
        read2 = createReadRecord(1, "1", 2800, 2900, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-300);

        geneSet.clearCounts();

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(FragmentType.UNSPLICED)]);

        // alternative splicing - first from reads with splits
        geneSet.clearCounts();

        read1 = createReadRecord(1, "1", 1050, 1100, REF_BASE_STR_1, createCigar(0, 50, 5000, 100, 0));
        read1.setFragmentInsertSize(500);
        read2 = createReadRecord(1, "1", 3100, 3300, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-500);

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(FragmentType.ALT)]);

        int longInsertSize = config.MaxFragmentLength + 100;

        // alt splicing - exon to exon read skipping an exon and long, currently not detected
        geneSet.clearCounts();

        read1 = createReadRecord(1, "1", 1050, 1100, REF_BASE_STR_1, createCigar(0, 50, 0));
        read1.setFragmentInsertSize(longInsertSize);
        read2 = createReadRecord(1, "1", 4550, 4650, REF_BASE_STR_1, createCigar(0, 100, 0));
        read2.setFragmentInsertSize(-longInsertSize);

        reads = Lists.newArrayList(read1, read2);
        bamReader.processReadRecords(geneSet, reads);

        assertEquals(1, geneCounts[typeAsInt(TOTAL)]);
        assertEquals(1, geneCounts[typeAsInt(FragmentType.ALT)]);
    }

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