package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.svtools.rna_expression.RnaBamReader.findStringOverlaps;
import static com.hartwig.hmftools.svtools.rna_expression.RnaBamReader.overlaps;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.junit.Test;

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
    public void testPositionMatches()
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
}