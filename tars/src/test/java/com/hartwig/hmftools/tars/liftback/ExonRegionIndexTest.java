package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createGeneDataCache;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.junit.Test;

// contains() uses binary search on merged ranges; earlier scan-back was O(N) and tripped on overlap edge cases.
public class ExonRegionIndexTest
{
    // builds an index from a single gene on ensembl chromosome "1" (keyed as chr1 under V38) whose one
    // transcript carries the given exon spans. fromCache merges the spans, so overlapping/abutting spans coalesce.
    private static ExonRegionIndex indexFor(final String chromosome, final List<int[]> exonSpans)
    {
        final EnsemblDataCache cache = createGeneDataCache();
        addGeneData(cache, chromosome, List.of(createEnsemblGeneData("ENSG_TEST", "TESTG", chromosome, 1, 1, 100_000)));

        final TranscriptData transcript = new TranscriptData(
                1, "ENST_TEST", "ENSG_TEST", true, (byte) 1, 1, 100_000, null, null, "protein_coding", "");
        final List<ExonData> exons = new ArrayList<>();
        int rank = 1;
        for(final int[] span : exonSpans)
            exons.add(new ExonData(1, span[0], span[1], rank++, -1, -1));
        transcript.setExons(exons);
        addTransExonData(cache, "ENSG_TEST", List.of(transcript));

        return ExonRegionIndex.fromCache(cache, V38);
    }

    @Test
    public void testContains()
    {
        final ExonRegionIndex idx = indexFor("1", List.of(new int[] { 100, 200 }, new int[] { 300, 400 }));

        // hits include interval boundaries
        assertTrue(idx.contains("chr1", 100));
        assertTrue(idx.contains("chr1", 150));
        assertTrue(idx.contains("chr1", 200));
        assertTrue(idx.contains("chr1", 300));
        assertTrue(idx.contains("chr1", 400));

        // intronic / intergenic positions miss
        assertFalse(idx.contains("chr1", 50));
        assertFalse(idx.contains("chr1", 250));
        assertFalse(idx.contains("chr1", 500));
    }

    @Test
    public void testContainsAcrossOverlappingAndAbuttingIntervals()
    {
        final ExonRegionIndex idx = indexFor("1", List.of(
                new int[] { 100, 200 },
                new int[] { 150, 300 },
                new int[] { 301, 400 },
                new int[] { 600, 700 }));

        assertTrue(idx.contains("chr1", 250));
        assertTrue(idx.contains("chr1", 301));
        assertTrue(idx.contains("chr1", 400));
        assertFalse(idx.contains("chr1", 401));
        assertFalse(idx.contains("chr1", 599));
        assertTrue(idx.contains("chr1", 700));
    }

    @Test
    public void testChromosomeKeyedInVersionedForm()
    {
        // ensembl stores bare "1"; the index keys it in the run's ref-genome form (V38 -> chr1) so lookups by
        // lifted genomic contig match directly. The bare form does not match.
        final ExonRegionIndex idx = indexFor("1", List.of(new int[] { 100, 200 }));
        assertTrue(idx.contains("chr1", 150));
        assertFalse(idx.contains("1", 150));
    }
}
