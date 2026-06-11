package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;

import org.junit.Test;

public class ContigTranslatorTest
{
    private static final String GENE_ID = "ENSG00000000001";
    private static final String GENE_NAME = "TESTGN";
    private static final String TRANS_NAME = "ENST00000000001";
    private static final String CONTIG_NAME = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;

    // three exons on chr1: 100-199, 300-399, 500-549; introns 200-299, 400-499.
    private static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                CONTIG_NAME, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    @Test
    public void testReadEntirelyInFirstExon()
    {
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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 201, cigar("30M"));

        assertNotNull(result);
        assertEquals(500, result.genomicStart());
        assertEquals("30M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadCrossingOneJunction()
    {
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
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("130M"));

        assertNotNull(result);
        assertEquals(190, result.genomicStart());
        assertEquals("10M100N100M100N20M", result.genomicCigar().toString());
        assertEquals(2, result.impliedIntrons().size());
        assertEquals(new BaseRegion(200, 299), result.impliedIntrons().get(0));
        assertEquals(new BaseRegion(400, 499), result.impliedIntrons().get(1));
    }

    @Test
    public void testReadWithSoftClip()
    {
        // 5S consumes query only — genomic start is pos-based, not query-based.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("5S100M"));

        assertNotNull(result);
        assertEquals(150, result.genomicStart());
        assertEquals("5S50M100N50M", result.genomicCigar().toString());
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
        // fills exon 1 exactly — guards against a spurious trailing N
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("100M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("100M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadExtendingPastLastSpanIsClampedToTrailingSoftClip()
    {
        // overhang past last span is converted to trailing soft-clip
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 231, cigar("50M"));

        assertNotNull(result);
        assertEquals(530, result.genomicStart());
        assertEquals("20M30S", result.genomicCigar().toString());
    }

    @Test
    public void testContigPosBeyondContigLengthReturnsNull()
    {
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
        assertNull(ContigTranslator.translate(threeExonContig(), -5, cigar("4M10S")));
    }

    @Test
    public void testReadAtFirstBaseOfContig()
    {
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("1M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("1M", result.genomicCigar().toString());
    }

    @Test
    public void testReadAtLastBaseOfContig()
    {
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 250, cigar("1M"));

        assertNotNull(result);
        assertEquals(549, result.genomicStart());
        assertEquals("1M", result.genomicCigar().toString());
    }

    @Test
    public void testReadWithLeadingAndTrailingSoftClips()
    {
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("10S100M10S"));

        assertNotNull(result);
        assertEquals(150, result.genomicStart());
        assertEquals("10S50M100N50M10S", result.genomicCigar().toString());
    }

    @Test
    public void testDeletionCrossingJunctionAboveThreshold()
    {
        // D above SPLICE_FLANKING_DELETION_MAX_BP is preserved after the intron N, not absorbed.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("10M10D5M"));

        assertNotNull(result);
        assertEquals(190, result.genomicStart());
        assertEquals("10M100N10D5M", result.genomicCigar().toString());
        assertEquals(1, result.impliedIntrons().size());
    }

    @Test
    public void testDeletionCrossingJunctionCollapsedIntoSplice()
    {
        // tx-FASTA off-by-N artefact: small D straddling an exon boundary is absorbed into the splice N.
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
        ContigEntry twoExon = new ContigEntry(
                "ensG_X_T", 1, 100, "G", "X", "T", CHR_1, 1,
                List.of(new BaseRegion(100, 149), new BaseRegion(200, 249)));

        // fills exon 1 exactly — no spurious trailing N
        ContigTranslator.TranslationResult result = ContigTranslator.translate(twoExon, 1, cigar("50M"));
        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("50M", result.genomicCigar().toString());

        result = ContigTranslator.translate(twoExon, 50, cigar("51M"));
        assertNotNull(result);
        assertEquals(149, result.genomicStart());
        assertEquals("1M50N50M", result.genomicCigar().toString());
    }

    @Test
    public void testSoftClipAtTrailingExonBoundary()
    {
        assertTrue(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 71, cigar("30M20S")));
    }

    @Test
    public void testSoftClipAtLeadingExonBoundary()
    {
        assertTrue(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 101, cigar("20S30M")));
    }

    @Test
    public void testSoftClipNotAtBoundary()
    {
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 60, cigar("30M20S")));
    }

    @Test
    public void testNoSoftClipReturnsFalse()
    {
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 71, cigar("30M")));
    }

    @Test
    public void testSingleExonContigReturnsFalse()
    {
        // no interior boundaries on a single-exon contig
        ContigEntry singleExon = new ContigEntry(
                "single", 1, 100, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1, 1, List.of(new BaseRegion(100, 199)));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(singleExon, 1, cigar("50S50M")));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(singleExon, 1, cigar("50M50S")));
    }

    @Test
    public void testContigOuterEdgesNotTreatedAsBoundary()
    {
        // outer edges of the contig are not interior boundaries
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 1, cigar("20S30M")));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 201, cigar("50M20S")));
    }

    @Test
    public void testTrimMicroAnchorAtThreshold()
    {
        final Cigar in = cigar("100M500N3M48S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 2);
        assertEquals("100M500N3M48S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorBelowThreshold()
    {
        final Cigar in = cigar("100M500N2M49S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("100M51S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorAtThresholdKept()
    {
        // 3M tail with threshold 3 is kept — trim fires only strictly below threshold.
        final Cigar in = cigar("100M500N3M48S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("100M500N3M48S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorOneBaseTail()
    {
        final Cigar in = cigar("80M2000N1M70S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("80M71S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorNoTrailingSoftClip()
    {
        final Cigar in = cigar("80M2000N2M");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("80M2000N2M", out.toString());
    }

    @Test
    public void testTrimMicroAnchorNoTrailingIntron()
    {
        // no preceding N -> nothing to trim
        final Cigar in = cigar("100M");
        assertEquals("100M", ContigTranslator.trimTrailingMicroAnchor(in, 3).toString());

        // adjacent M ops are not a junction
        final Cigar in2 = cigar("80M3M48S");
        assertEquals("80M3M48S", ContigTranslator.trimTrailingMicroAnchor(in2, 3).toString());
    }

    @Test
    public void testTrimMicroAnchorMustHaveAnchorBeforeIntron()
    {
        // no preceding M anchor — trimming would leave only softclips, so refuse.
        final Cigar in = cigar("50S100N2M99S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("50S100N2M99S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorTrailingMatchAboveThresholdUnchanged()
    {
        final Cigar in = cigar("100M500N4M47S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("100M500N4M47S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorsLeading()
    {
        // 2M anchor dropped, rolled into leading S; start advances by anchor + intron = 3715.
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("38S2M3713N68M4I39M"), 1, 3);
        assertEquals("40S68M4I39M", out.AdjustedCigar.toString());
        assertEquals(3715, out.StartShift);
    }

    @Test
    public void testTrimMicroAnchorsLeadingAnchorAtThresholdKept()
    {
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("38S3M3713N68M"), 1, 3);
        assertEquals("38S3M3713N68M", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);
    }

    @Test
    public void testTrimMicroAnchorsBothEnds()
    {
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("10S2M500N80M600N1M20S"), 1, 3);
        assertEquals("12S80M21S", out.AdjustedCigar.toString());
        assertEquals(502, out.StartShift);
    }

    @Test
    public void testTrimMicroAnchorsNoChange()
    {
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("38S68M3713N40M"), 1, 3);
        assertEquals("38S68M3713N40M", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);
    }

    @Test
    public void testTrimMicroAnchorsBareTinyAnchorKept()
    {
        // bare leading anchor (no softclip) means the read genuinely starts mid-exon — keep it.
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("2M80N149M"), 1, 3);
        assertEquals("2M80N149M", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);
    }

    @Test
    public void testTrimMicroAnchorsSoftclipAdjacentTinyAnchorDropped()
    {
        // same 2bp anchor but with trailing softclip: bwa over-ran the exon boundary, junction unsupported — drop it.
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("117M89N2M32S"), 1, 3);
        assertEquals("117M34S", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);
    }

    private static Cigar cigar(final String cigarString)
    {
        return TextCigarCodec.decode(cigarString);
    }
}
