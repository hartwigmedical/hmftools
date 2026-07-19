package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.exonRegionIndex;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

// contains() uses binary search on merged ranges; earlier scan-back was O(N) and tripped on overlap edge cases.
public class ExonRegionIndexTest
{
    @Test
    public void testContains()
    {
        ExonRegionIndex exonIndex = exonRegionIndex("1", List.of(new int[] { 100, 200 }, new int[] { 300, 400 }));

        // hits include interval boundaries
        assertTrue(exonIndex.contains("chr1", 100));
        assertTrue(exonIndex.contains("chr1", 150));
        assertTrue(exonIndex.contains("chr1", 200));
        assertTrue(exonIndex.contains("chr1", 300));
        assertTrue(exonIndex.contains("chr1", 400));

        // intronic / intergenic positions miss
        assertFalse(exonIndex.contains("chr1", 50));
        assertFalse(exonIndex.contains("chr1", 250));
        assertFalse(exonIndex.contains("chr1", 500));
    }

    @Test
    public void testContainsAcrossOverlappingAndAbuttingIntervals()
    {
        ExonRegionIndex exonIndex = exonRegionIndex("1", List.of(
                new int[] { 100, 200 },
                new int[] { 150, 300 },
                new int[] { 301, 400 },
                new int[] { 600, 700 }));

        assertTrue(exonIndex.contains("chr1", 250));
        assertTrue(exonIndex.contains("chr1", 301));
        assertTrue(exonIndex.contains("chr1", 400));
        assertFalse(exonIndex.contains("chr1", 401));
        assertFalse(exonIndex.contains("chr1", 599));
        assertTrue(exonIndex.contains("chr1", 700));
    }

    @Test
    public void testChromosomeKeyedInVersionedForm()
    {
        // ensembl stores bare "1"; the index keys it in the run's ref-genome form (V38 -> chr1) so lookups by
        // lifted genomic contig match directly. The bare form does not match.
        ExonRegionIndex exonIndex = exonRegionIndex("1", List.of(new int[] { 100, 200 }));
        assertTrue(exonIndex.contains("chr1", 150));
        assertFalse(exonIndex.contains("1", 150));
    }
}
