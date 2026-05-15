package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
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
                CONTIG_NAME, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    @Test
    public void testReadEntirelyInFirstExon()
    {
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("50M"));

        assertNotNull(result);
        assertEquals(CHR_1, result.Chromosome);
        assertEquals(100, result.GenomicStart);
        assertEquals("50M", result.GenomicCigar.toString());
        assertTrue(result.ImpliedIntrons.isEmpty());
    }

    @Test
    public void testReadEntirelyInLastExon()
    {
        // contig pos 201 is within the third span (cumulative offset is 200 after spans 1+2)
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 201, cigar("30M"));

        assertNotNull(result);
        assertEquals(500, result.GenomicStart);
        assertEquals("30M", result.GenomicCigar.toString());
        assertTrue(result.ImpliedIntrons.isEmpty());
    }

    @Test
    public void testReadCrossingOneJunction()
    {
        // start at contig pos 51 (genomic 150), 100M extends 100 bases. First 50 are in exon 1 (150-199),
        // remaining 50 fall into exon 2 starting at 300. Intron is 200-299 (length 100).
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("100M"));

        assertNotNull(result);
        assertEquals(150, result.GenomicStart);
        assertEquals("50M100N50M", result.GenomicCigar.toString());
        assertEquals(1, result.ImpliedIntrons.size());
        assertEquals(new BaseRegion(200, 299), result.ImpliedIntrons.get(0));
    }

    @Test
    public void testReadCrossingTwoJunctions()
    {
        // start at contig pos 91 (genomic 190), 130M = 10 in exon1 + 100 in exon2 + 20 in exon3
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("130M"));

        assertNotNull(result);
        assertEquals(190, result.GenomicStart);
        assertEquals("10M100N100M100N20M", result.GenomicCigar.toString());
        assertEquals(2, result.ImpliedIntrons.size());
        assertEquals(new BaseRegion(200, 299), result.ImpliedIntrons.get(0));
        assertEquals(new BaseRegion(400, 499), result.ImpliedIntrons.get(1));
    }

    @Test
    public void testReadWithSoftClip()
    {
        // 5S100M starting at contig pos 51. The 5S consumes query only — genomic start should still be 150.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("5S100M"));

        assertNotNull(result);
        assertEquals(150, result.GenomicStart);
        assertEquals("5S50M100N50M", result.GenomicCigar.toString());
    }

    @Test
    public void testReadWithInsertion()
    {
        // 30M5I20M starting at contig pos 1 — entirely in exon 1, no junctions
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("30M5I20M"));

        assertNotNull(result);
        assertEquals(100, result.GenomicStart);
        assertEquals("30M5I20M", result.GenomicCigar.toString());
        assertTrue(result.ImpliedIntrons.isEmpty());
    }

    @Test
    public void testReadWithDeletion()
    {
        // 30M5D20M starting at contig pos 1 — D consumes ref but doesn't cross a junction (ends at genomic 154)
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("30M5D20M"));

        assertNotNull(result);
        assertEquals(100, result.GenomicStart);
        assertEquals("30M5D20M", result.GenomicCigar.toString());
        assertTrue(result.ImpliedIntrons.isEmpty());
    }

    @Test
    public void testReadExactlyFillingFirstExon()
    {
        // 100M starting at contig pos 1 — fills exon 1 exactly (100-199), no spurious N at the end
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("100M"));

        assertNotNull(result);
        assertEquals(100, result.GenomicStart);
        assertEquals("100M", result.GenomicCigar.toString());
        assertTrue(result.ImpliedIntrons.isEmpty());
    }

    @Test
    public void testReadExtendingPastLastSpanIsUnliftable()
    {
        // start at contig pos 231 (genomic 530), 50M extends past the last span (500-549, only 20 bases remain)
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 231, cigar("50M"));

        assertNull(result);
    }

    @Test
    public void testContigPosBeyondContigLengthReturnsNull()
    {
        // contig length is 250 (100+100+50). Pos 251 is past the end.
        assertNull(ContigTranslator.translate(threeExonContig(), 251, cigar("10M")));
    }

    @Test
    public void testInvalidContigPosReturnsNull()
    {
        assertNull(ContigTranslator.translate(threeExonContig(), 0, cigar("10M")));
        assertNull(ContigTranslator.translate(threeExonContig(), -5, cigar("10M")));
    }

    @Test
    public void testReadAtFirstBaseOfContig()
    {
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("1M"));

        assertNotNull(result);
        assertEquals(100, result.GenomicStart);
        assertEquals("1M", result.GenomicCigar.toString());
    }

    @Test
    public void testReadAtLastBaseOfContig()
    {
        // contig length is 250, so pos 250 + 1M lands on the very last contig base
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 250, cigar("1M"));

        assertNotNull(result);
        assertEquals(549, result.GenomicStart);
        assertEquals("1M", result.GenomicCigar.toString());
    }

    @Test
    public void testReadWithLeadingAndTrailingSoftClips()
    {
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("10S100M10S"));

        assertNotNull(result);
        assertEquals(150, result.GenomicStart);
        assertEquals("10S50M100N50M10S", result.GenomicCigar.toString());
    }

    @Test
    public void testDeletionCrossingJunction()
    {
        // 10M5D5M starting at contig pos 91 (genomic 190): 10M fills the rest of exon 1, then 5D would consume
        // ref starting at 200 (intron). We emit 100N for the intron, advance to exon 2 (300), then continue 5D + 5M.
        ContigTranslator.TranslationResult result = ContigTranslator.translate(threeExonContig(), 91, cigar("10M5D5M"));

        assertNotNull(result);
        assertEquals(190, result.GenomicStart);
        assertEquals("10M100N5D5M", result.GenomicCigar.toString());
        assertEquals(1, result.ImpliedIntrons.size());
    }

    @Test
    public void testTwoExonContigBoundaryCases()
    {
        ContigEntry twoExon = new ContigEntry(
                "ensG_X_T", "G", "X", "T", CHR_1,
                List.of(new BaseRegion(100, 149), new BaseRegion(200, 249)));

        // M that exactly fills exon 1 — no spurious trailing N
        ContigTranslator.TranslationResult result = ContigTranslator.translate(twoExon, 1, cigar("50M"));
        assertNotNull(result);
        assertEquals(100, result.GenomicStart);
        assertEquals("50M", result.GenomicCigar.toString());

        // last base of exon 1 + 51M reaches one base into exon 2
        result = ContigTranslator.translate(twoExon, 50, cigar("51M"));
        assertNotNull(result);
        assertEquals(149, result.GenomicStart);
        assertEquals("1M50N50M", result.GenomicCigar.toString());
    }

    private static Cigar cigar(final String cigarString)
    {
        return TextCigarCodec.decode(cigarString);
    }
}
