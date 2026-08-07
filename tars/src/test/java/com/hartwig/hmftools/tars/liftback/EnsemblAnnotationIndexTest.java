package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

public class EnsemblAnnotationIndexTest
{
    private static final String CHR_1 = "chr1";
    private static final String CHR_2 = "chr2";

    @Test
    public void testContainsExon()
    {
        EnsemblAnnotationIndex annotationIndex = annotationIndex("1", List.of(new int[] { 100, 200 }, new int[] { 300, 400 }));

        assertTrue(annotationIndex.containsExon("chr1", 100));
        assertTrue(annotationIndex.containsExon("chr1", 150));
        assertTrue(annotationIndex.containsExon("chr1", 200));
        assertTrue(annotationIndex.containsExon("chr1", 300));
        assertTrue(annotationIndex.containsExon("chr1", 400));

        assertFalse(annotationIndex.containsExon("chr1", 50));
        assertFalse(annotationIndex.containsExon("chr1", 250));
        assertFalse(annotationIndex.containsExon("chr1", 500));
    }

    @Test
    public void testContainsAcrossOverlappingAndAbuttingIntervals()
    {
        EnsemblAnnotationIndex annotationIndex = annotationIndex(
                "1", List.of(
                        new int[] { 100, 200 },
                        new int[] { 150, 300 },
                        new int[] { 301, 400 },
                        new int[] { 600, 700 }));

        assertTrue(annotationIndex.containsExon("chr1", 250));
        assertTrue(annotationIndex.containsExon("chr1", 301));
        assertTrue(annotationIndex.containsExon("chr1", 400));
        assertFalse(annotationIndex.containsExon("chr1", 401));
        assertFalse(annotationIndex.containsExon("chr1", 599));
        assertTrue(annotationIndex.containsExon("chr1", 700));
    }

    @Test
    public void testNestedExonKeepsEnclosingSpan()
    {
        // an exon nested inside an earlier-starting one must not shrink the merged span: the union end wins, so
        // the enclosing exon's tail (past the nested exon's end) still reports contained.
        EnsemblAnnotationIndex annotationIndex = annotationIndex("1", List.of(new int[] { 100, 200 }, new int[] { 100, 150 }));
        assertTrue(annotationIndex.containsExon("chr1", 175));
        assertTrue(annotationIndex.containsExon("chr1", 200));
        assertFalse(annotationIndex.containsExon("chr1", 201));
    }

    @Test
    public void testChromosomeKeyedInVersionedForm()
    {
        // ensembl stores bare "1"; the index keys chromosomes in the run's ref-genome form (V38 -> chr1).
        EnsemblAnnotationIndex annotationIndex = annotationIndex("1", List.of(new int[] { 100, 200 }));
        assertTrue(annotationIndex.containsExon("chr1", 150));
        assertFalse(annotationIndex.containsExon("1", 150));
    }

    @Test
    public void testContainsJunctionByValue()
    {
        EnsemblAnnotationIndex annotationIndex = junctionIndex(intron(CHR_1, 200, 299));

        assertTrue(annotationIndex.containsJunction(intron(CHR_1, 200, 299)));
        assertFalse(annotationIndex.containsJunction(intron(CHR_1, 200, 300)));
        assertFalse(annotationIndex.containsJunction(intron(CHR_2, 200, 299)));
    }

    @Test
    public void testJunctionBoundaryLookups()
    {
        ChrBaseRegion shortIntron = intron(CHR_1, 200, 299);
        ChrBaseRegion longIntron = intron(CHR_1, 200, 450);
        ChrBaseRegion farIntron = intron(CHR_1, 50, 299);
        EnsemblAnnotationIndex annotationIndex = junctionIndex(shortIntron, longIntron, farIntron);

        assertTrue(annotationIndex.junctionsByStart(CHR_1, 200).contains(shortIntron));
        assertTrue(annotationIndex.junctionsByStart(CHR_1, 200).contains(longIntron));
        assertTrue(annotationIndex.junctionsByEnd(CHR_1, 299).contains(shortIntron));
        assertTrue(annotationIndex.junctionsByEnd(CHR_1, 299).contains(farIntron));
        assertTrue(annotationIndex.junctionsByStart(CHR_2, 200).isEmpty());
        assertTrue(annotationIndex.junctionsByEnd(CHR_1, 500).isEmpty());
    }

    @Test
    public void testEmptyJunctionIndex()
    {
        EnsemblAnnotationIndex annotationIndex = EnsemblAnnotationIndex.fromJunctions((Set<ChrBaseRegion>) null);

        assertEquals(0, annotationIndex.junctionCount());
        assertFalse(annotationIndex.containsJunction(intron(CHR_1, 200, 299)));
        assertTrue(annotationIndex.junctionsByStart(CHR_1, 200).isEmpty());
    }

    @Test
    public void testJunctionsFromAdjacentExons()
    {
        EnsemblAnnotationIndex annotationIndex = EnsemblAnnotationIndex.fromContigEntries(List.of(
                entry("chr1", new int[] { 1000, 1099 }, new int[] { 1500, 1599 }, new int[] { 2000, 2099 }),
                entry("chr1", new int[] { 1000, 1099 }, new int[] { 1700, 1799 })));

        assertEquals(3, annotationIndex.junctionCount());
        assertTrue(annotationIndex.containsJunction(new ChrBaseRegion("chr1", 1100, 1499)));
        assertTrue(annotationIndex.containsJunction(new ChrBaseRegion("chr1", 1600, 1999)));
        assertTrue(annotationIndex.containsJunction(new ChrBaseRegion("chr1", 1100, 1699)));
        assertEquals(1, annotationIndex.junctionStrand(new ChrBaseRegion("chr1", 1100, 1499)));
    }

    @Test
    public void testJunctionSharedByOppositeStrandsIsUnknown()
    {
        ContigEntry forward = entry("chr2", 1, new int[] { 100, 199 }, new int[] { 300, 399 });
        ContigEntry reverse = entry("chr2", -1, new int[] { 100, 199 }, new int[] { 300, 399 });

        EnsemblAnnotationIndex annotationIndex = EnsemblAnnotationIndex.fromContigEntries(List.of(forward, reverse));

        assertEquals(0, annotationIndex.junctionStrand(new ChrBaseRegion("chr2", 200, 299)));
    }

    @Test
    public void testSingleExonEntryProducesNoJunctions()
    {
        EnsemblAnnotationIndex annotationIndex = EnsemblAnnotationIndex.fromContigEntries(List.of(
                entry("chr1", new int[] { 1000, 1999 })));

        assertEquals(0, annotationIndex.junctionCount());
    }

    @Test
    public void testDuplicateJunctionsAcrossEntriesDeduplicate()
    {
        EnsemblAnnotationIndex annotationIndex = EnsemblAnnotationIndex.fromContigEntries(List.of(
                entry("chr2", new int[] { 100, 199 }, new int[] { 300, 399 }),
                entry("chr2", new int[] { 100, 199 }, new int[] { 300, 399 })));

        assertEquals(1, annotationIndex.junctionCount());
        assertTrue(annotationIndex.containsJunction(new ChrBaseRegion("chr2", 200, 299)));
    }

    @Test
    public void testExonsOutOfOrderStillGiveLowToHighJunction()
    {
        EnsemblAnnotationIndex annotationIndex = EnsemblAnnotationIndex.fromContigEntries(List.of(
                entry("chr4", new int[] { 2000, 2099 }, new int[] { 1000, 1099 })));

        assertEquals(1, annotationIndex.junctionCount());
        assertTrue(annotationIndex.containsJunction(new ChrBaseRegion("chr4", 1100, 1999)));
    }

    private static EnsemblAnnotationIndex annotationIndex(final String chromosome, final List<int[]> exonSpans)
    {
        List<BaseRegion> spans = new ArrayList<>();
        for(int[] span : exonSpans)
            spans.add(new BaseRegion(span[0], span[1]));

        ContigEntry entry = ContigEntry.annotationOnly("gene", "GENE", "trans", V38.versionedChromosome(chromosome), 1, spans);
        return EnsemblAnnotationIndex.fromContigEntries(List.of(entry));
    }

    private static EnsemblAnnotationIndex junctionIndex(final ChrBaseRegion... introns)
    {
        return EnsemblAnnotationIndex.fromJunctions(new HashSet<>(List.of(introns)));
    }

    private static ChrBaseRegion intron(final String chromosome, final int start, final int end)
    {
        return new ChrBaseRegion(chromosome, start, end);
    }

    private static ContigEntry entry(final String chromosome, final int[]... exons)
    {
        return entry(chromosome, 1, exons);
    }

    private static ContigEntry entry(final String chromosome, final int strand, final int[]... exons)
    {
        List<BaseRegion> spans = new ArrayList<>();
        for(int[] exon : exons)
        {
            spans.add(new BaseRegion(exon[0], exon[1]));
        }
        return ContigEntry.annotationOnly("g", "gn", "tn", chromosome, strand, spans);
    }
}
