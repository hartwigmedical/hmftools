package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

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
        ExonRegionIndex exonIndex = exonRegionIndex(
                "1", List.of(
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
    public void testNestedExonKeepsEnclosingSpan()
    {
        // an exon nested inside an earlier-starting one must not shrink the merged span: the union end wins, so
        // the enclosing exon's tail (past the nested exon's end) still reports contained.
        ExonRegionIndex exonIndex = exonRegionIndex("1", List.of(new int[] { 100, 200 }, new int[] { 100, 150 }));
        assertTrue(exonIndex.contains("chr1", 175));
        assertTrue(exonIndex.contains("chr1", 200));
        assertFalse(exonIndex.contains("chr1", 201));
    }

    @Test
    public void testChromosomeKeyedInVersionedForm()
    {
        // ensembl stores bare "1"; the index keys chromosomes in the run's ref-genome form (V38 -> chr1).
        ExonRegionIndex exonIndex = exonRegionIndex("1", List.of(new int[] { 100, 200 }));
        assertTrue(exonIndex.contains("chr1", 150));
        assertFalse(exonIndex.contains("1", 150));
    }

    private static ExonRegionIndex exonRegionIndex(final String chromosome, final List<int[]> exonSpans)
    {
        List<BaseRegion> spans = new ArrayList<>();
        for(int[] span : exonSpans)
            spans.add(new BaseRegion(span[0], span[1]));

        ContigEntry entry = ContigEntry.annotationOnly("gene", "GENE", "trans", V38.versionedChromosome(chromosome), 1, spans);
        return ExonRegionIndex.fromContigEntries(List.of(entry));
    }
}
