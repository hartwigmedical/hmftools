package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

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

// All tests use TarsTestFixtures.threeExonContig(), whose packed contig coordinates map to the genome as:
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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("50M"));

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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 201, cigar("30M"));

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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("100M"));

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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("130M"));

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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("30M5I20M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("30M5I20M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadWithDeletion()
    {
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("30M5D20M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("30M5D20M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadExactlyFillingFirstExon()
    {
        // fills exon 1 exactly - guards against a spurious trailing N
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("100M"));

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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 231, cigar("50M"));

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
        // pos 0 is 1 base before altStart; leading M absorbs into soft-clip.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 0, cigar("10M"));

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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 250, cigar("1M"));

        assertNotNull(result);
        assertEquals(549, result.genomicStart());
        assertEquals("1M", result.genomicCigar().toString());
    }

    @Test
    public void testReadWithLeadingAndTrailingSoftClips()
    {
        // softclips consume query only - genomic start is pos-based, not query-based; both clips survive the lift.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("10S100M10S"));

        assertNotNull(result);
        assertEquals(150, result.genomicStart());
        assertEquals("10S50M100N50M10S", result.genomicCigar().toString());
    }

    @Test
    public void testDeletionCrossingJunctionAboveThreshold()
    {
        // contig 91 -> chr1 190; the 10M reaches the exon1 end (contig 100), so the D sits right at the
        // boundary. 10D exceeds SPLICE_FLANKING_DELETION_MAX_BP (5) so it is preserved after the intron N, not absorbed.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("10M10D5M"));

        assertNotNull(result);
        assertEquals(190, result.genomicStart());
        assertEquals("10M100N10D5M", result.genomicCigar().toString());
        assertEquals(1, result.impliedIntrons().size());
    }

    @Test
    public void testDeletionCrossingJunctionCollapsedIntoSplice()
    {
        // tx-FASTA off-by-N artefact: 5D (at SPLICE_FLANKING_DELETION_MAX_BP = 5) straddling an exon boundary is absorbed into the splice N.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("10M5D5M"));

        assertNotNull(result);
        assertEquals(190, result.genomicStart());
        assertEquals("10M105N5M", result.genomicCigar().toString());
        assertEquals(1, result.impliedIntrons().size());
    }

    @Test
    public void testDeletionsBothSidesOfJunctionCollapsedIntoSplice()
    {
        // D straddles the boundary: lift splits it, then both flanking Ds are absorbed into the N.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 92, cigar("9M3D6M"));

        assertNotNull(result);
        assertEquals(191, result.genomicStart());
        assertEquals("9M103N6M", result.genomicCigar().toString());
        assertEquals(1, result.impliedIntrons().size());
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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(twoExon, 1, cigar("50M"));
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
    public void testHasSoftClipAtExonBoundary()
    {
        // softclip abutting an interior exon boundary, trailing then leading -> true
        assertTrue(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 71, cigar("30M20S")));
        assertTrue(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 101, cigar("20S30M")));

        // softclip away from any boundary, and no softclip at all -> false
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 60, cigar("30M20S")));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 71, cigar("30M")));

        // single-exon contig has no interior boundaries -> false either side
        ContigEntry singleExon = new ContigEntry(
                "single", 1, 100, "G", "X", "T", CHR_1, 1, List.of(new BaseRegion(100, 199)));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(singleExon, 1, cigar("50S50M")));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(singleExon, 1, cigar("50M50S")));

        // outer edges of the contig are not interior boundaries -> false
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 1, cigar("20S30M")));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 201, cigar("50M20S")));
    }

    @Test
    public void testTrimTrailingMicroAnchor()
    {
        // 3M tail at threshold 2 -> kept (trim fires only strictly below threshold)
        assertEquals("100M500N3M48S", ContigTranslator.trimTrailingMicroAnchor(cigar("100M500N3M48S"), 2).toString());

        // 2M tail below threshold 3 -> dropped, intron collapsed into the trailing clip
        assertEquals("100M51S", ContigTranslator.trimTrailingMicroAnchor(cigar("100M500N2M49S"), 3).toString());

        // 3M tail at threshold 3 -> kept
        assertEquals("100M500N3M48S", ContigTranslator.trimTrailingMicroAnchor(cigar("100M500N3M48S"), 3).toString());

        // 1bp tail -> dropped
        assertEquals("80M71S", ContigTranslator.trimTrailingMicroAnchor(cigar("80M2000N1M70S"), 3).toString());

        // no trailing softclip to fold into -> unchanged
        assertEquals("80M2000N2M", ContigTranslator.trimTrailingMicroAnchor(cigar("80M2000N2M"), 3).toString());

        // no preceding N -> nothing to trim; adjacent M ops are not a junction
        assertEquals("100M", ContigTranslator.trimTrailingMicroAnchor(cigar("100M"), 3).toString());
        assertEquals("80M3M48S", ContigTranslator.trimTrailingMicroAnchor(cigar("80M3M48S"), 3).toString());

        // no preceding M anchor - trimming would leave only softclips, so refuse
        assertEquals("50S100N2M99S", ContigTranslator.trimTrailingMicroAnchor(cigar("50S100N2M99S"), 3).toString());

        // 4M tail above threshold -> kept
        assertEquals("100M500N4M47S", ContigTranslator.trimTrailingMicroAnchor(cigar("100M500N4M47S"), 3).toString());
    }

    @Test
    public void testTrimMicroAnchors()
    {
        // leading 2M anchor dropped, rolled into leading S; start advances by anchor + intron = 3715
        ContigTranslator.MicroAnchorResult out = ContigTranslator.trimMicroAnchors(cigar("38S2M3713N68M4I39M"), 1, 3);
        assertEquals("40S68M4I39M", out.AdjustedCigar.toString());
        assertEquals(3715, out.StartShift);

        // leading 3M anchor at threshold -> kept
        out = ContigTranslator.trimMicroAnchors(cigar("38S3M3713N68M"), 1, 3);
        assertEquals("38S3M3713N68M", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);

        // both ends trimmed
        out = ContigTranslator.trimMicroAnchors(cigar("10S2M500N80M600N1M20S"), 1, 3);
        assertEquals("12S80M21S", out.AdjustedCigar.toString());
        assertEquals(502, out.StartShift);

        // nothing tiny -> unchanged
        out = ContigTranslator.trimMicroAnchors(cigar("38S68M3713N40M"), 1, 3);
        assertEquals("38S68M3713N40M", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);

        // bare leading anchor (no softclip) means the read genuinely starts mid-exon -> keep it
        out = ContigTranslator.trimMicroAnchors(cigar("2M80N149M"), 1, 3);
        assertEquals("2M80N149M", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);

        // same 2bp anchor but with trailing softclip: bwa over-ran the boundary, junction unsupported -> drop it
        out = ContigTranslator.trimMicroAnchors(cigar("117M89N2M32S"), 1, 3);
        assertEquals("117M34S", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);
    }

    @Test
    public void testContiguousSpansEmitNoZeroLengthIntron()
    {
        // adjacent spans (no gap) imply no intron, so no 0-length N element is emitted
        ContigEntry contiguous = new ContigEntry(
                "ensG_X_T", 1, 100, "G", "X", "T", CHR_1, 1,
                List.of(new BaseRegion(100, 149), new BaseRegion(150, 199)));

        ContigTranslator.TranslationResult result = ContigTranslator.translate(contiguous, 1, cigar("100M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertFalse(result.genomicCigar().getCigarElements().stream().anyMatch(x -> x.getLength() == 0));
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

    private static Cigar cigar(final String cigarString)
    {
        return TextCigarCodec.decode(cigarString);
    }
}
