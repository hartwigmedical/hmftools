package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;

// threeExonContig(), used by most tests below, maps packed contig coordinates to the genome as:
//
//   contig 1..100   -> exon1  chr1 100..199   (100 bp)
//   contig 101..200 -> exon2  chr1 300..399   (100 bp)
//   contig 201..250 -> exon3  chr1 500..549   ( 50 bp)
//   introns (implied as N on lift): 200..299 (100 bp) and 400..499 (100 bp)
//
// So a contig position lifts to genomic = exon.start + (contigPos - exon.contigStart). E.g. contig 51 -> 150
// (exon1), contig 91 -> 190 (exon1), contig 201 -> 500 (exon3), contig 250 -> 549 (last base).
public class ContigTranslatorTest
{
    @Test
    public void testReadEntirelyInFirstExon()
    {
        // contig 1..50 sits inside exon1 -> chr1 100, no boundary crossed.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("50M"));

        assertNotNull(result);
        assertEquals(CHR_1, result.chromosome());
        assertEquals(100, result.genomicStart());
        assertEquals("50M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadEntirelyInLastExon()
    {
        // contig 201..230 sits inside exon3 -> chr1 500, no boundary crossed.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 201, cigar("30M"));

        assertNotNull(result);
        assertEquals(500, result.genomicStart());
        assertEquals("30M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadCrossingOneJunction()
    {
        // contig 51..150 spans the exon1/exon2 boundary (contig 100|101): 50 bp of exon1 (from chr1 150),
        // then intron1 as 100N, then 50 bp of exon2.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("100M"));

        assertNotNull(result);
        assertEquals(150, result.genomicStart());
        assertEquals("50M100N50M", result.genomicCigar().toString());
        assertEquals(1, result.impliedIntrons().size());
        assertEquals(new BaseRegion(200, 299), result.impliedIntrons().get(0));
    }

    @Test
    public void testReadCrossingTwoJunctions()
    {
        // contig 91..220 spans both boundaries: 10 bp tail of exon1 (from chr1 190), intron1 (100N), all
        // 100 bp of exon2, intron2 (100N), 20 bp head of exon3.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("130M"));

        assertNotNull(result);
        assertEquals(190, result.genomicStart());
        assertEquals("10M100N100M100N20M", result.genomicCigar().toString());
        assertEquals(2, result.impliedIntrons().size());
        assertEquals(new BaseRegion(200, 299), result.impliedIntrons().get(0));
        assertEquals(new BaseRegion(400, 499), result.impliedIntrons().get(1));
    }

    @Test
    public void testReadWithInsertion()
    {
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("30M5I20M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("30M5I20M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadWithDeletion()
    {
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("30M5D20M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("30M5D20M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadExactlyFillingFirstExon()
    {
        // fills exon 1 exactly - guards against a spurious trailing N
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("100M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("100M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadExtendingPastLastSpanIsClampedToTrailingSoftClip()
    {
        // contig 231 -> exon3 chr1 530; only 20 bp remain to the contig end (250), so 20M then the 30 bp
        // overhang past the last span becomes a trailing soft-clip.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 231, cigar("50M"));

        assertNotNull(result);
        assertEquals(530, result.genomicStart());
        assertEquals("20M30S", result.genomicCigar().toString());
    }

    @Test
    public void testContigPosBeyondContigLengthReturnsNull()
    {
        // contig 251 is past the contig end (length 250) -> untranslatable.
        assertNull(ContigTranslator.translate(threeExonContig(), 251, cigar("10M")));
    }

    @Test
    public void testLeadingOverhangClampedToSoftClip()
    {
        // pos 0 is 1 base before contigStart; leading M absorbs into soft-clip.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 0, cigar("10M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("1S9M", result.genomicCigar().toString());
    }

    @Test
    public void testLeadingOverhangExceedsLeadingMReturnsNull()
    {
        // contig start -5 is 6 bases before the contig (pos 1), but only 4M lead it -> overhang can't be
        // absorbed into a soft-clip, so the lift fails.
        assertNull(ContigTranslator.translate(threeExonContig(), -5, cigar("4M10S")));
    }

    @Test
    public void testReadAtLastBaseOfContig()
    {
        // contig 250 is the final base of exon3 -> chr1 549.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 250, cigar("1M"));

        assertNotNull(result);
        assertEquals(549, result.genomicStart());
        assertEquals("1M", result.genomicCigar().toString());
    }

    @Test
    public void testReadWithLeadingAndTrailingSoftClips()
    {
        // leading softclip does not shift genomic start (pos-based, not query-based); both clips survive the lift.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("10S100M10S"));

        assertNotNull(result);
        assertEquals(150, result.genomicStart());
        assertEquals("10S50M100N50M10S", result.genomicCigar().toString());
    }

    @Test
    public void testDeletionCrossingBoundaryLeftIntact()
    {
        // a deletion abutting the exon1/exon2 junction is emitted intact after the intron N; translate neither splits nor absorbs it.
        ContigTranslateResult small = ContigTranslator.translate(threeExonContig(), 91, cigar("10M5D5M"));
        assertNotNull(small);
        assertEquals(190, small.genomicStart());
        assertEquals("10M100N5D5M", small.genomicCigar().toString());
        assertEquals(1, small.impliedIntrons().size());

        ContigTranslateResult large = ContigTranslator.translate(threeExonContig(), 91, cigar("10M10D5M"));
        assertEquals("10M100N10D5M", large.genomicCigar().toString());
    }

    @Test
    public void testTwoExonContigBoundaryCases()
    {
        // a different two-exon contig: contig 1..50 -> exon1 chr1 100..149, contig 51..100 -> exon2 chr1 200..249;
        // implied intron 150..199.
        ContigEntry twoExon = new ContigEntry(
                "ensG_X_T", 1, 100, "G", "X", "T", CHR_1, 1,
                List.of(new BaseRegion(100, 149), new BaseRegion(200, 249)));

        // contig 1..50 fills exon1 exactly -> chr1 100, no spurious trailing N.
        ContigTranslateResult result = ContigTranslator.translate(twoExon, 1, cigar("50M"));
        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("50M", result.genomicCigar().toString());

        // contig 50..100 crosses the boundary: 1 bp tail of exon1 (chr1 149), 50N intron, 50 bp of exon2.
        result = ContigTranslator.translate(twoExon, 50, cigar("51M"));
        assertNotNull(result);
        assertEquals(149, result.genomicStart());
        assertEquals("1M50N50M", result.genomicCigar().toString());
    }

    @Test
    public void testSoftClipAtExonBoundaryReported()
    {
        // softclip abutting an interior exon boundary, trailing then leading -> true
        assertTrue(ContigTranslator.translate(threeExonContig(), 71, cigar("30M20S")).softClipAtExonBoundary());
        assertTrue(ContigTranslator.translate(threeExonContig(), 101, cigar("20S30M")).softClipAtExonBoundary());

        // softclip away from any boundary, and no softclip at all -> false
        assertFalse(ContigTranslator.translate(threeExonContig(), 60, cigar("30M20S")).softClipAtExonBoundary());
        assertFalse(ContigTranslator.translate(threeExonContig(), 71, cigar("30M")).softClipAtExonBoundary());

        // single-exon contig has no interior boundaries -> false either side
        ContigEntry singleExon = new ContigEntry(
                "single", 1, 100, "G", "X", "T", CHR_1, 1, List.of(new BaseRegion(100, 199)));
        assertFalse(ContigTranslator.translate(singleExon, 1, cigar("50S50M")).softClipAtExonBoundary());
        assertFalse(ContigTranslator.translate(singleExon, 1, cigar("50M50S")).softClipAtExonBoundary());

        // outer edges of the contig are not interior boundaries -> false
        assertFalse(ContigTranslator.translate(threeExonContig(), 1, cigar("20S30M")).softClipAtExonBoundary());
        assertFalse(ContigTranslator.translate(threeExonContig(), 201, cigar("50M20S")).softClipAtExonBoundary());
    }

    @Test
    public void testContiguousSpansEmitNoZeroLengthIntron()
    {
        // adjacent spans (no gap) imply no intron, so no 0-length N element is emitted
        ContigEntry contiguous = new ContigEntry(
                "ensG_X_T", 1, 100, "G", "X", "T", CHR_1, 1,
                List.of(new BaseRegion(100, 149), new BaseRegion(150, 199)));

        ContigTranslateResult result = ContigTranslator.translate(contiguous, 1, cigar("100M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertFalse(result.genomicCigar().getCigarElements().stream().anyMatch(element -> element.getLength() == 0));
        assertEquals("100M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testDropZeroLengthElementsAndMergeFlankingIntrons()
    {
        // a 0-length interior element (eg a 0M for a zero-span exon between two introns) is dropped, then the
        // flanking introns merge into one
        List<CigarElement> elements = cigar("72M102N0M2538N79M").getCigarElements();
        Cigar result = new Cigar(ContigTranslator.mergeAdjacentSameOp(ContigTranslator.dropZeroLength(elements)));
        assertEquals("72M2640N79M", result.toString());
    }

    private static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                "ensG_T_tx", 1, 250, "G", "GENE", "T", CHR_1, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    private static Cigar cigar(final String cigarString)
    {
        return TextCigarCodec.decode(cigarString);
    }
}
