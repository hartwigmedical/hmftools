package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.common.TarsCigarUtils;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;

// threeExonContig(), used by most tests below, maps packed contig coordinates to the genome as:
//
//   contig 1..100   -> exon1  chr1 100..199   (100 bp)
//   contig 101..200 -> exon2  chr1 300..399   (100 bp)
//   contig 201..250 -> exon3  chr1 500..549   ( 50 bp)
//   introns (implied as N on lift): 200..299 (100 bp) and 400..499 (100 bp)
//
// A contig position lifts to genomic = exon.start + (contigPos - exon.contigStart).
public class ContigTranslatorTest
{
    @Test
    public void testReadEntirelyInFirstExon()
    {
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
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 201, cigar("30M"));

        assertNotNull(result);
        assertEquals(500, result.genomicStart());
        assertEquals("30M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadCrossingOneJunction()
    {
        // contig 51..150 spans the exon1/exon2 boundary: 50 bp of exon1 from chr1 150, then intron1 as 100N, then 50 bp of exon2.
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
        // contig 91..220 spans both boundaries: 10 bp of exon1 from chr1 190, intron1 100N, all of exon2, intron2 100N, 20 bp of exon3.
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
        // guards against a spurious trailing N
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 1, cigar("100M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("100M", result.genomicCigar().toString());
        assertTrue(result.impliedIntrons().isEmpty());
    }

    @Test
    public void testReadExtendingPastLastSpanIsClampedToTrailingSoftClip()
    {
        // only 20 bp remain to contigEnd(250), so the 30 bp overhang past the last span becomes a trailing soft-clip.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 231, cigar("50M"));

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
        // pos 0 is one base before contigStart.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 0, cigar("10M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("1S9M", result.genomicCigar().toString());
    }

    @Test
    public void testLeadingOverhangExceedsLeadingMReturnsNull()
    {
        // 6 bases before contigStart with only 4M leading, so the overhang cannot be absorbed into a soft-clip.
        assertNull(ContigTranslator.translate(threeExonContig(), -5, cigar("4M10S")));
    }

    @Test
    public void testLeadingOverhangAbsorbedAcrossDeletion()
    {
        // 11 bases before contigStart with only 7M leading: the 4D supplies the remaining 4 reference bases, so the
        // clamp still lands on contigStart. Seen on real data as 7M4D144M.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), -10, cigar("7M4D80M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("7S80M", result.genomicCigar().toString());
    }

    @Test
    public void testLeadingOverhangAbsorbedAcrossInsertion()
    {
        // 8 bases before contigStart, led by 6M1I: the insertion adds read bases but no reference, so 2 more come from
        // the next M and all three elements fold into the clip. Seen on real data as 6M1I144M.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), -7, cigar("6M1I80M"));

        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("9S78M", result.genomicCigar().toString());
    }

    @Test
    public void testTrailingOverhangAbsorbedAcrossDeletion()
    {
        // read ends 5 bases past contigEnd(250) with only 4M trailing; the preceding 1D supplies the fifth. Seen on
        // real data as 147M1D4M.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 211, cigar("40M1D4M"));

        assertNotNull(result);
        assertEquals(510, result.genomicStart());
        assertEquals("40M4S", result.genomicCigar().toString());
    }

    @Test
    public void testLeadingOverhangEndingOnDeletionReturnsNull()
    {
        // the overhang is absorbed exactly by 3M, leaving a deletion at the clip boundary. Dropping it would start the
        // alignment past contigStart, so the lift is declined rather than misplaced.
        assertNull(ContigTranslator.translate(threeExonContig(), -2, cigar("3M1D80M")));
    }

    @Test
    public void testReadAtLastBaseOfContig()
    {
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 250, cigar("1M"));

        assertNotNull(result);
        assertEquals(549, result.genomicStart());
        assertEquals("1M", result.genomicCigar().toString());
    }

    @Test
    public void testReadWithLeadingAndTrailingSoftClips()
    {
        // a leading soft-clip does not shift the genomic start; both clips survive the lift.
        ContigTranslateResult result = ContigTranslator.translate(threeExonContig(), 51, cigar("10S100M10S"));

        assertNotNull(result);
        assertEquals(150, result.genomicStart());
        assertEquals("10S50M100N50M10S", result.genomicCigar().toString());
    }

    @Test
    public void testDeletionCrossingBoundaryLeftIntact()
    {
        // translate() emits a deletion abutting the junction intact after the N; absorbing it is mergeDeletionsIntoSplice's job.
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
        // contig 1..50 -> chr1 100..149, contig 51..100 -> chr1 200..249; implied intron 150..199.
        ContigEntry twoExon = new ContigEntry(
                "ensG_X_T", 1, 100, "G", "X", "T", CHR_1, 1,
                List.of(new BaseRegion(100, 149), new BaseRegion(200, 249)));

        // fills exon1 exactly: no spurious trailing N
        ContigTranslateResult result = ContigTranslator.translate(twoExon, 1, cigar("50M"));
        assertNotNull(result);
        assertEquals(100, result.genomicStart());
        assertEquals("50M", result.genomicCigar().toString());

        result = ContigTranslator.translate(twoExon, 50, cigar("51M"));
        assertNotNull(result);
        assertEquals(149, result.genomicStart());
        assertEquals("1M50N50M", result.genomicCigar().toString());
    }

    @Test
    public void testXaAlignmentsAreLiftedAndDeduplicated()
    {
        ContigTranslator translator = new ContigTranslator(List.of(threeExonContig()));

        List<LiftedAlignment> alignments = translator.liftXaAlignments(
                "ensG_T_tx,+51,50M,1;ensG_T_tx,+51,50M,1;ensG_T_tx,-201,30M,2;bad;unknown_tx,+1,10M,0;");

        assertEquals(2, alignments.size());
        assertEquals(150, alignments.get(0).LiftedPos);
        assertEquals("50M", alignments.get(0).LiftedCigar);
        assertTrue(alignments.get(0).ForwardStrand);
        assertEquals(500, alignments.get(1).LiftedPos);
        assertFalse(alignments.get(1).ForwardStrand);
    }

    @Test
    public void testContiguousSpansEmitNoZeroLengthIntron()
    {
        // adjacent spans imply no intron, so no 0-length N element is emitted
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
        // a 0-length interior element (a 0M for a zero-span exon between two introns) is dropped, then the flanking introns merge
        List<CigarElement> elements = cigar("72M102N0M2538N79M").getCigarElements();
        Cigar result = new Cigar(TarsCigarUtils.normalize(elements));
        assertEquals("72M2640N79M", result.toString());
    }

    @Test
    public void testAbsorbsSmallDeletionAfterJunction()
    {
        // 5D is within MAX_MERGED_DELETION_BP, so it folds into the N
        assertEquals("10M105N5M", mergeDeletions("10M100N5D5M"));
    }

    @Test
    public void testAbsorbsSmallDeletionBeforeJunction()
    {
        assertEquals("10M105N5M", mergeDeletions("10M5D100N5M"));
    }

    @Test
    public void testAbsorbsDeletionsBothSidesOfJunction()
    {
        assertEquals("10M105N5M", mergeDeletions("10M3D100N2D5M"));
    }

    @Test
    public void testKeepsDeletionAboveThreshold()
    {
        // 6D exceeds MAX_MERGED_DELETION_BP: a real deletion, kept beside the N
        assertEquals("10M100N6D5M", mergeDeletions("10M100N6D5M"));
    }

    @Test
    public void testLeavesPlainSpliceUntouched()
    {
        assertEquals("50M100N50M", mergeDeletions("50M100N50M"));
    }

    @Test
    public void testDoesNotFoldALeadingDeletionIntoTheSplice()
    {
        // with no aligned block ahead of it the D must not be absorbed: a cigar cannot begin with N, and 9760N149M
        // reached a written BAM exactly this way
        assertEquals("2D9758N149M", mergeDeletions("2D9758N149M"));
    }

    @Test
    public void testDoesNotFoldATrailingDeletionIntoTheSplice()
    {
        // mirror of the leading case: absorbing the final D would end the cigar on an N, which isofox reads as a splice
        // junction with no exon after it
        assertEquals("149M9758N2D", mergeDeletions("149M9758N2D"));
    }

    private static String mergeDeletions(final String cigarString)
    {
        return new Cigar(ContigTranslator.mergeDeletionsIntoSplice(cigar(cigarString).getCigarElements())).toString();
    }

    private static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                "ensG_T_tx", 1, 250, "G", "GENE", "T", CHR_1, 1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

    private static Cigar cigar(final String cigarString)
    {
        return CigarUtils.cigarFromStr(cigarString);
    }
}
