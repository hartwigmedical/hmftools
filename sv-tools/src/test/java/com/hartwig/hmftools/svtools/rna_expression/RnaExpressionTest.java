package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.svtools.rna_expression.RnaBamReader.findStringOverlaps;
import static com.hartwig.hmftools.svtools.rna_expression.RnaBamReader.overlaps;
import static com.hartwig.hmftools.svtools.rna_expression.RnaBamReader.setMatchingBases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

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

        ReadRecord record = new ReadRecord("1", "1", 800, 900, "", null);
        assertFalse(overlaps(region, record));

        record = new ReadRecord("1", "1", 2200, 2300, "", null);
        assertFalse(overlaps(region, record));

        record = new ReadRecord("1", "1", 1500, 1600, "", null);
        assertFalse(overlaps(region, record));

        record = new ReadRecord("1", "1", 800, 2200, "", null);
        assertFalse(overlaps(region, record));

        record = new ReadRecord("1", "1", 800, 1200, "", null);
        assertTrue(overlaps(region, record));

        record = new ReadRecord("1", "1", 1900, 2200, "", null);
        assertTrue(overlaps(region, record));
    }

    @Test
    public void testMatchingBases()
    {
        String refBaseString = "ABCDEFGHIJKLMNOPQRST";
        String refBaseString2 = "TSRQPONMLKJIHGFEDCBA";

        // test a read covering part of a region
        RegionReadData region = new RegionReadData(GenomeRegions.create("1", 100, 119), "1");
        region.setRefBases(refBaseString);

        ReadRecord record = new ReadRecord("1", "1", 105, 114, refBaseString.substring(5, 15), createCigar(0, 10, 0));

        setMatchingBases(region, record);

        assertEquals(10, region.baseCoverage(1));
        assertEquals(1, region.matchedReadCount());

        // a read covering all of an exon with unmapped bases either end
        region = new RegionReadData(GenomeRegions.create("1", 105, 114), "1");
        region.setRefBases(refBaseString.substring(5, 15));

        record = new ReadRecord("1", "1", 105, 124, refBaseString, createCigar(5, 10, 5));
        setMatchingBases(region, record);

        assertEquals(10, region.baseCoverage(1));
        assertEquals(1, region.matchedReadCount());

        // soft-clipped ends with matching position at one end
        record = new ReadRecord("1", "1", 105, 119, refBaseString.substring(2, 19),
                createCigar(3, 10, 5, 0, 0));

        region.clearState();
        setMatchingBases(region, record);

        assertEquals(10, region.baseCoverage(1));
        assertEquals(1, region.matchedReadCount());

        // read overlapping 2 different regions
        region = new RegionReadData(GenomeRegions.create("1", 100, 119), "1");
        region.setRefBases(refBaseString);

        // first match 105-114 = 10 bases, then 100 non-reference, then 115-214, then 215-224
        record = new ReadRecord("1", "1", 105, 224, refBaseString.substring(5, 15) + refBaseString2.substring(0, 10),
                createCigar(0, 10, 100, 10, 0));

        setMatchingBases(region, record);

        assertEquals(10, region.baseCoverage(1));
        assertEquals(1, region.matchedReadCount());

        // matching exact exon boundary at start
        region = new RegionReadData(GenomeRegions.create("1", 105, 114), "1");
        region.setRefBases(refBaseString.substring(5, 15));

        setMatchingBases(region, record);

        assertEquals(10, region.baseCoverage(1));
        assertEquals(1, region.matchedReadCount());

        // and a region at the other end
        region = new RegionReadData(GenomeRegions.create("1", 215, 224), "1");
        region.setRefBases(refBaseString2.substring(0, 10));

        setMatchingBases(region, record);

        assertEquals(10, region.baseCoverage(1));
        assertEquals(1, region.matchedReadCount());

        // a region extending further than the read at the other end, with the first 5 bases matching a region before the exon
        String intronicBases = "ABCDE";
        record = new ReadRecord("1", "1", 105, 224, refBaseString.substring(5, 15) + intronicBases + refBaseString2.substring(0, 5),
                createCigar(0, 10, 100, 10, 0));

        region = new RegionReadData(GenomeRegions.create("1", 220, 229), "1");
        region.setRefBases(refBaseString2.substring(0, 10));

        setMatchingBases(region, record);

        assertEquals(5, region.baseCoverage(1));
        assertEquals(0, region.matchedReadCount()); // insufficient to count as a match
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