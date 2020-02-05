package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_MATCH;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.WITHIN_EXON;
import static com.hartwig.hmftools.svtools.rna_expression.RnaBamReader.findStringOverlaps;
import static com.hartwig.hmftools.svtools.rna_expression.RnaBamReader.overlaps;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.ALT;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.EXONIC;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.UNSPLICED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class RnaExpressionTest
{
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
    public void testPositionOverlaps()
    {
        GenomeRegion region = GenomeRegions.create("1", 1000, 2000);

        ReadRecord record = createReadRecord(1, "1", 800, 900, "", null);
        assertFalse(overlaps(region, record));

        record = createReadRecord(1, "1", 2200, 2300, "", null);
        assertFalse(overlaps(region, record));

        record = createReadRecord(1, "1", 1500, 1600, "", null);
        assertFalse(overlaps(region, record));

        record = createReadRecord(1, "1", 800, 2200, "", null);
        assertFalse(overlaps(region, record));

        record = createReadRecord(1, "1", 800, 1200, "", null);
        assertTrue(overlaps(region, record));

        record = createReadRecord(1, "1", 1900, 2200, "", null);
        assertTrue(overlaps(region, record));
    }

    @Test
    public void testReadRegionTypes()
    {
        String refBaseString = "ABCDEFGHIJKLMNOPQRST"; // currently unused

        // single read and exon
        RegionReadData region = new RegionReadData(GenomeRegions.create("1", 100, 200));

        // read covers part of the exon
        ReadRecord read = createReadRecord(1, "1", 110, 130, refBaseString, createCigar(0, 21, 0));

        assertEquals(1, read.getMappedRegionCoords().size());
        assertEquals(110, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(130, read.getMappedRegionCoords().get(0)[SE_END]);

        // test classification of reads
        assertEquals(WITHIN_EXON, read.getRegionMatchType(region));

        read = createReadRecord(1, "1", 90, 110, refBaseString, createCigar(0, 21, 0));
        assertEquals(EXON_INTRON, read.getRegionMatchType(region));

        read = createReadRecord(1, "1", 100, 200, refBaseString, createCigar(0, 101, 0));
        assertEquals(EXON_MATCH, read.getRegionMatchType(region));

        read = createReadRecord(1, "1", 100, 150, refBaseString, createCigar(0, 51, 0));
        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region));

        // read covering multiple exons
        Cigar cigar = new Cigar();
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches half of first exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(21, CigarOperator.M)); // matches all of second exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(25, CigarOperator.M)); // matches past last exon into intron

        read = createReadRecord(1, "1", 110, 204, refBaseString, cigar);

        assertEquals(3, read.getMappedRegionCoords().size());
        assertEquals(110, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(120, read.getMappedRegionCoords().get(0)[SE_END]);
        assertEquals(140, read.getMappedRegionCoords().get(1)[SE_START]);
        assertEquals(160, read.getMappedRegionCoords().get(1)[SE_END]);
        assertEquals(180, read.getMappedRegionCoords().get(2)[SE_START]);
        assertEquals(204, read.getMappedRegionCoords().get(2)[SE_END]);

        RegionReadData region1 = new RegionReadData(GenomeRegions.create("1", 100, 120));
        RegionReadData region2 = new RegionReadData(GenomeRegions.create("1", 140, 160));
        RegionReadData region3 = new RegionReadData(GenomeRegions.create("1", 180, 200));

        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region1));
        assertEquals(EXON_MATCH, read.getRegionMatchType(region2));
        assertEquals(EXON_INTRON, read.getRegionMatchType(region3));
    }

    @Test
    public void testUnmappedSpliceJunctions()
    {
        String refBaseString = "ABCDEFGHIJKLMNOPQRST"; // currently unused

        // soft-clipping at one end
        RegionReadData region1 = new RegionReadData(GenomeRegions.create("1", 141, 150));
        region1.setRefBases(refBaseString.substring(0, 10));
        RegionReadData region2 = new RegionReadData(GenomeRegions.create("1", 200, 209));
        region2.setRefBases(refBaseString.substring(10, 19));

        region1.addPostRegion(region2);
        region2.addPreRegion(region1);

        ReadRecord read = createReadRecord(1, "1", 200, 209, refBaseString.substring(5, 19), createCigar(5, 10, 0));

        List<RegionReadData> regions = Lists.newArrayList(region2);
        read.processOverlappingRegions(regions);

        assertEquals(EXON_BOUNDARY, read.getRegionMatchType(region1));
        assertEquals(2, read.getMappedRegionCoords().size());
        assertEquals(146, read.getMappedRegionCoords().get(0)[SE_START]);
        assertEquals(150, read.getMappedRegionCoords().get(0)[SE_END]);

        // test again on the up side, with 2 different regions matching the inferred bases
        read = createReadRecord(1, "1", 141, 155, refBaseString.substring(0, 15), createCigar(0, 15, 0));

        RegionReadData region3 = new RegionReadData(GenomeRegions.create("1", 200, 229));
        region3.setRefBases(refBaseString.substring(10, 19) + refBaseString);

        RegionReadData region4 = new RegionReadData(GenomeRegions.create("1", 121, 150));
        region1.setRefBases(refBaseString + refBaseString.substring(0, 10));

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
    }

    @Test
    public void testReadTranscriptClassification()
    {
        String refBaseString = "ABCDEFGHIJKLMNOPQRST"; // currently unused

        String trans1 = "TRANS01";

        RegionReadData region = createRegion(trans1, 1, "1", 100, 200);

        // unspliced read does not support the transcript
        ReadRecord read = createReadRecord(1, "1", 90, 110, refBaseString, createCigar(0, 21, 0));

        List<RegionReadData> regions = Lists.newArrayList(region);
        read.processOverlappingRegions(regions);

        assertEquals(UNSPLICED, read.getTranscriptClassification(trans1));

        // exonic read does support the transcript
        read = createReadRecord(1, "1", 120, 140, refBaseString, createCigar(0, 21, 0));

        regions = Lists.newArrayList(region);
        read.processOverlappingRegions(regions);

        assertEquals(EXONIC, read.getTranscriptClassification(trans1));

        // skipped exon
        RegionReadData region1 = createRegion(trans1,1,"1", 100, 120);
        RegionReadData region2 = createRegion(trans1,2, "1", 140, 160);
        RegionReadData region3 = createRegion(trans1,3, "1", 180, 200);

        // read covering multiple exons but skips the middle exon
        read = createReadRecord(1, "1", 110, 200, refBaseString, createCigar(0, 11, 59, 21, 0));

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

        read = createReadRecord(1, "1", 110, 200, refBaseString, cigar);

        read.processOverlappingRegions(regions);

        assertEquals(ALT, read.getTranscriptClassification(trans1));

        // read spanning into intro (unspliced)
        cigar = new Cigar();
        cigar.add(new CigarElement(11, CigarOperator.M)); // matches half of first exon
        cigar.add(new CigarElement(19, CigarOperator.N));
        cigar.add(new CigarElement(31, CigarOperator.M)); // goes past second exon into intron
        cigar.add(new CigarElement(9, CigarOperator.N));
        cigar.add(new CigarElement(25, CigarOperator.M)); // matches past last exon into intron

        read = createReadRecord(1, "1", 110, 200, refBaseString, cigar);

        read.processOverlappingRegions(regions);

        assertEquals(UNSPLICED, read.getTranscriptClassification(trans1));
    }

    private ReadRecord createReadRecord(
            final int id, final String chromosome, long posStart, long posEnd, final String readBases, final Cigar cigar)
    {
        Cigar readCigar = cigar != null ? cigar : createCigar(0, (int)(posEnd - posStart + 1), 0);
        return new ReadRecord(null, String.valueOf(id), chromosome, posStart, posEnd, readBases, readCigar);
    }

    private RegionReadData createRegion(final String trans, int exonRank, final String chromosome, long posStart, long posEnd)
    {
        RegionReadData region = new RegionReadData(GenomeRegions.create(chromosome, posStart, posEnd));
        region.addExonRef(trans, exonRank);
        return region;
    }

    @Test
    public void testCigarCreation()
    {
        Cigar cigar = createCigar(2, 10,1);
        assertTrue(cigar.toString().equals("2S10M1S"));

        cigar = createCigar(0, 10,100, 12, 0);
        assertTrue(cigar.toString().equals("10M100N12M"));

        cigar = createCigar(2, 10,100, 12, 4);
        assertTrue(cigar.toString().equals("2S10M100N12M4S"));
    }

    public static Cigar createCigar(int preSc, int match, int postSc)
    {
        if(preSc == 0 && match == 0 &&postSc == 0)
            return null;

        Cigar cigar = new Cigar();

        if(preSc > 0)
            cigar.add(new CigarElement(preSc, CigarOperator.SOFT_CLIP));

        if(match > 0)
            cigar.add(new CigarElement(match, CigarOperator.MATCH_OR_MISMATCH));

        if(postSc > 0)
            cigar.add(new CigarElement(postSc, CigarOperator.SOFT_CLIP));

        return cigar;

    }

    public static Cigar createCigar(int preSc, int preMatch, int nonRef, int postMatch, int postSc)
    {
        if(preSc == 0 && preMatch == 0 && nonRef == 0 && postMatch == 0 && postSc == 0)
            return null;

        Cigar cigar = new Cigar();

        if(preSc > 0)
            cigar.add(new CigarElement(preSc, CigarOperator.SOFT_CLIP));

        if(preMatch > 0)
            cigar.add(new CigarElement(preMatch, CigarOperator.MATCH_OR_MISMATCH));

        if(nonRef > 0)
            cigar.add(new CigarElement(nonRef, CigarOperator.SKIPPED_REGION));

        if(postMatch > 0)
            cigar.add(new CigarElement(postMatch, CigarOperator.MATCH_OR_MISMATCH));

        if(postSc > 0)
            cigar.add(new CigarElement(postSc, CigarOperator.SOFT_CLIP));

        return cigar;
    }

}