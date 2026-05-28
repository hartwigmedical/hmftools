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

    // three exon spans on chr1: 100-199 (len 100), 300-399 (len 100), 500-549 (len 50). Two introns: 200-299, 400-499.
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
        // contig pos 201 is within the third span (cumulative offset is 200 after spans 1+2)
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 201, cigar("30M"));

        assertNotNull(result);
        assertEquals(500, result.genomicStart());
        assertEquals("30M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadCrossingOneJunction()
    {
        // start at contig pos 51 (genomic 150), 100M extends 100 bases. First 50 are in exon 1 (150-199),
        // remaining 50 fall into exon 2 starting at 300. Intron is 200-299 (length 100).
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
        // start at contig pos 91 (genomic 190), 130M = 10 in exon1 + 100 in exon2 + 20 in exon3
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
        // 5S100M starting at contig pos 51. The 5S consumes query only — genomic start should still be 150.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("5S100M"));

        assertNotNull(result);
        assertEquals(150, result.genomicStart());
        assertEquals("5S50M100N50M", result.genomicCigar().toString());
    }

    @Test
    public void testReadWithInsertion()
    {
        // 30M5I20M starting at contig pos 1 — entirely in exon 1, no junctions
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("30M5I20M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("30M5I20M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadWithDeletion()
    {
        // 30M5D20M starting at contig pos 1 — D consumes ref but doesn't cross a junction (ends at genomic 154)
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("30M5D20M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("30M5D20M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadExactlyFillingFirstExon()
    {
        // 100M starting at contig pos 1 — fills exon 1 exactly (100-199), no spurious N at the end
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("100M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("100M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadExtendingPastLastSpanIsClampedToTrailingSoftClip()
    {
        // start at contig pos 231 (genomic 530), 50M extends past the last span (500-549, only 20 bases remain).
        // overhang of 30 bases is converted into trailing soft-clip.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 231, cigar("50M"));

        assertNotNull(result);
        assertEquals(530, result.genomicStart());
        assertEquals("20M30S", result.genomicCigar().toString());
    }

    @Test
    public void testContigPosBeyondContigLengthReturnsNull()
    {
        // contig length is 250 (100+100+50). Pos 251 is past the end — no overhang fix possible.
        assertNull(ContigTranslator.translate(threeExonContig(), 251, cigar("10M")));
    }

    @Test
    public void testLeadingOverhangClampedToSoftClip()
    {
        // contig pos 0 is 1 base before altStart (1); leading M absorbs the shift into soft-clip.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 0, cigar("10M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("1S9M", result.genomicCigar().toString());
    }

    @Test
    public void testLeadingOverhangExceedsLeadingMReturnsNull()
    {
        // pos -5 is 6 below altStart, but only 4M to absorb -> unliftable
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
        // contig length is 250, so pos 250 + 1M lands on the very last contig base
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
        // 10M10D5M starting at contig pos 91 (genomic 190): 10M fills the rest of exon 1, then 10D would
        // consume ref starting at 200 (intron). We emit 100N for the intron, advance to exon 2 (300),
        // then continue 10D + 5M. 10D is above SPLICE_FLANKING_DELETION_MAX_BP so it is preserved as-is.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("10M10D5M"));

        assertNotNull(result);
        assertEquals(190, result.genomicStart());
        assertEquals("10M100N10D5M", result.genomicCigar().toString());
        assertEquals(1, result.impliedIntrons().size());
    }

    @Test
    public void testDeletionCrossingJunctionCollapsedIntoSplice()
    {
        // tx-FASTA off-by-N artefact: bwa reports e.g. 10M5D5M where the 5D straddles an exon boundary.
        // The lift would emit 10M 100N 5D 5M; collapseSpliceFlankingDeletions absorbs the small trailing
        // D into the N because each D <= SPLICE_FLANKING_DELETION_MAX_BP.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("10M5D5M"));

        assertNotNull(result);
        assertEquals(190, result.genomicStart());
        assertEquals("10M105N5M", result.genomicCigar().toString());
        assertEquals(1, result.impliedIntrons().size());
    }

    @Test
    public void testDeletionsBothSidesOfJunctionCollapsedIntoSplice()
    {
        // pre-lift cigar with a single D straddling the boundary: 9M3D6M starting at contig pos 92.
        // 9M reaches the end of exon 1, the 3D splits into 1D before the junction + 2D after (during the
        // lift walk), and collapseSpliceFlankingDeletions then absorbs both into the N.
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

        // M that exactly fills exon 1 — no spurious trailing N
        ContigTranslator.TranslationResult result = ContigTranslator.translate(twoExon, 1, cigar("50M"));
        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("50M", result.genomicCigar().toString());

        // last base of exon 1 + 51M reaches one base into exon 2
        result = ContigTranslator.translate(twoExon, 50, cigar("51M"));
        assertNotNull(result);
        assertEquals(149, result.genomicStart());
        assertEquals("1M50N50M", result.genomicCigar().toString());
    }

    @Test
    public void testSoftClipAtTrailingExonBoundary()
    {
        // pos 71 + 30M ends at contig pos 100 = end of exon 1 (interior boundary). 20S abuts that boundary.
        assertTrue(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 71, cigar("30M20S")));
    }

    @Test
    public void testSoftClipAtLeadingExonBoundary()
    {
        // pos 101 = first base of exon 2 (cum length of exon 1 is 100). Leading 20S abuts that boundary.
        assertTrue(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 101, cigar("20S30M")));
    }

    @Test
    public void testSoftClipNotAtBoundary()
    {
        // 30M ends mid-exon-1 (contig pos 89), trailing S not at any interior boundary
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
        // single-exon contig has no interior boundaries; should always be false even with S
        ContigEntry singleExon = new ContigEntry(
                "single", 1, 100, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1, 1, List.of(new BaseRegion(100, 199)));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(singleExon, 1, cigar("50S50M")));
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(singleExon, 1, cigar("50M50S")));
    }

    @Test
    public void testContigOuterEdgesNotTreatedAsBoundary()
    {
        // pos 1 (start of first exon) + leading S — outer edge, not an interior boundary
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 1, cigar("20S30M")));
        // 50M from pos 201 ends at contig pos 250 = last base of last exon — outer edge, not interior
        assertFalse(ContigTranslator.hasSoftClipAtExonBoundary(threeExonContig(), 201, cigar("50M20S")));
    }

    @Test
    public void testTrimMicroAnchorAtThreshold()
    {
        // "100M500N3M48S" with min anchor 3 -> 3M stays (3 is not below 3).
        final Cigar in = cigar("100M500N3M48S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 2);
        assertEquals("100M500N3M48S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorBelowThreshold()
    {
        // "100M500N2M49S" with threshold 3 -> drop the N+2M, roll into trailing S: "100M51S".
        final Cigar in = cigar("100M500N2M49S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("100M51S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorAtThresholdKept()
    {
        // 3M tail with min anchor 3 is kept (STAR keeps anchors >= alignSJDBoverhangMin; trim < 3).
        final Cigar in = cigar("100M500N3M48S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("100M500N3M48S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorOneBaseTail()
    {
        // tiniest possible tail anchor (1M) gets rolled.
        final Cigar in = cigar("80M2000N1M70S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("80M71S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorNoTrailingSoftClip()
    {
        // No trailing S -> nothing to roll, unchanged.
        final Cigar in = cigar("80M2000N2M");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("80M2000N2M", out.toString());
    }

    @Test
    public void testTrimMicroAnchorNoTrailingIntron()
    {
        // Tail is yM zS but no preceding N -> nothing to drop, leave alone.
        final Cigar in = cigar("100M");
        assertEquals("100M", ContigTranslator.trimTrailingMicroAnchor(in, 3).toString());

        final Cigar in2 = cigar("80M3M48S");
        // adjacent M ops aren't a junction; trim must not fire.
        assertEquals("80M3M48S", ContigTranslator.trimTrailingMicroAnchor(in2, 3).toString());
    }

    @Test
    public void testTrimMicroAnchorMustHaveAnchorBeforeIntron()
    {
        // pathological case: cigar starts with S then N then tiny-M then S — no preceding M anchor,
        // trimming would leave just softclips. Refuse.
        final Cigar in = cigar("50S100N2M99S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("50S100N2M99S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorTrailingMatchAboveThresholdUnchanged()
    {
        // 4M tail with threshold 3 -> not trimmed (above threshold).
        final Cigar in = cigar("100M500N4M47S");
        final Cigar out = ContigTranslator.trimTrailingMicroAnchor(in, 3);
        assertEquals("100M500N4M47S", out.toString());
    }

    @Test
    public void testTrimMicroAnchorsLeading()
    {
        // "38S2M3713N68M4I39M" with min anchor 3 -> drop the leading 2M + 3713N, roll 2M into the
        // leading softclip: "40S68M4I39M". Start advances by 2 (anchor) + 3713 (intron) = 3715.
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("38S2M3713N68M4I39M"), 3);
        assertEquals("40S68M4I39M", out.AdjustedCigar.toString());
        assertEquals(3715, out.StartShift);
    }

    @Test
    public void testTrimMicroAnchorsLeadingAnchorAtThresholdKept()
    {
        // 3M leading anchor is kept (>= 3); no start shift.
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("38S3M3713N68M"), 3);
        assertEquals("38S3M3713N68M", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);
    }

    @Test
    public void testTrimMicroAnchorsBothEnds()
    {
        // tiny anchors at both ends: leading 2M and trailing 1M both dropped.
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("10S2M500N80M600N1M20S"), 3);
        assertEquals("12S80M21S", out.AdjustedCigar.toString());
        assertEquals(502, out.StartShift);
    }

    @Test
    public void testTrimMicroAnchorsNoChange()
    {
        // clean cigar with adequate anchors -> unchanged, no shift.
        final ContigTranslator.MicroAnchorResult out =
                ContigTranslator.trimMicroAnchors(cigar("38S68M3713N40M"), 3);
        assertEquals("38S68M3713N40M", out.AdjustedCigar.toString());
        assertEquals(0, out.StartShift);
    }

    private static Cigar cigar(final String cigarString)
    {
        return TextCigarCodec.decode(cigarString);
    }
}
