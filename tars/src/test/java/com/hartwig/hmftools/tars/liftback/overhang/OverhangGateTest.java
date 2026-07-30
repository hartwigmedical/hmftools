package com.hartwig.hmftools.tars.liftback.overhang;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

// Unit tests for the overhang gate: the keep condition, the trigger, the iterative junction collapse and the standalone softclip extension.
public class OverhangGateTest
{
    private static OverhangGate gate(final TestGenome genome)
    {
        return new OverhangGate(genome.asRefGenome());
    }

    private static TestGenome genome()
    {
        return new TestGenome().with(CHR_1, 5000, 'A');
    }

    @Test
    public void testKeepsAnchorScoringAboveFive()
    {
        // trailing "20M 100N 6M 4S": the 6bp anchor matches its far exon (score +6 > 5), so the junction is kept.
        // near exon chr1:1-20, intron 21-120, far exon anchor at chr1:121-126.
        TestGenome genome = genome().set(CHR_1, 121, "CCCCCC");
        OverhangGate.Result result = gate(genome).collapseJunctions(CHR_1, 1, "20M100N6M4S", bases("C".repeat(30)));

        assertEquals("20M100N6M4S", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testCollapsesAnchorScoringExactlyFive()
    {
        // trailing "20M 100N 5M 5S": the 5bp anchor matches its far exon (score +5), which does NOT clear the
        // strict > 5 bar, so the junction collapses; the intron read-through is all 'A' vs 'C' so nothing extends.
        TestGenome genome = genome().set(CHR_1, 121, "CCCCC");
        OverhangGate.Result result = gate(genome).collapseJunctions(CHR_1, 1, "20M100N5M5S", bases("C".repeat(30)));

        assertEquals("20M10S", result.cigar());
        assertTrue(result.dropped());
    }

    @Test
    public void testCollapsesShortAnchorOnScore()
    {
        // trailing "20M 100N 8M 4S": an 8bp anchor is still short (<= MIN_OVERHANG_LENGTH of 12) so it is scored - it
        // mismatches its far exon (all 'A' vs 'C') so it collapses. Length alone does not save a short anchor.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "20M100N8M4S", bases("C".repeat(32)));

        assertEquals("20M12S", result.cigar());
        assertTrue(result.dropped());
    }

    @Test
    public void testTrustsLongOverhangWithoutScoring()
    {
        // trailing "20M 100N 15M 4S": a 15bp anchor is longer than MIN_OVERHANG_LENGTH (12), so it is trusted outright
        // and the junction is kept even though those bases mismatch the far exon (all 'A' vs 'C', which would score < 5).
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "20M100N15M4S", bases("C".repeat(39)));

        assertEquals("20M100N15M4S", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testLeavesSingleJunctionNoSoftclipUntouched()
    {
        // "20M 100N 4M": a single junction with no adjacent softclip does not trigger the gate, so the 4bp anchor is never scored.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "20M100N4M", bases("C".repeat(24)));

        assertEquals("20M100N4M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testCollapsesMultiJunctionWhenReferenceScoresBetter()
    {
        // "20M 50N 5M 60N 4M" (no soft clip, two junctions trigger the gate). The read is all 'A'; the terminal 4M's
        // far exon (chr1:136-139) is set to 'C' so the overhang scores -16, while reading contiguously through the 60N
        // (chr1:76-79) matches (+4). The intronic reference scores higher than the overhang, so the 60N collapses: the
        // 4M reads into the 5M -> "20M 50N 9M". The remaining single junction has no soft clip, so the collapse stops.
        TestGenome genome = genome().set(CHR_1, 136, "CCCC");
        OverhangGate.Result result = gate(genome).collapseJunctions(CHR_1, 1, "20M50N5M60N4M", bases("A".repeat(29)));

        assertEquals("20M50N9M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testKeepsMultiJunctionWhenOverhangScoresPositive()
    {
        // "20M 50N 5M 60N 4M", read and genome all 'A': the terminal 4M aligns positively across its junction (+4 > 0),
        // so the junction is kept without even checking the reference alternative.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "20M50N5M60N4M", bases("A".repeat(29)));

        assertEquals("20M50N5M60N4M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testKeepsMultiJunctionWhenReferenceNotBetter()
    {
        // "20M 50N 5M 60N 4M", read all 'C' vs an all-'A' genome: the terminal 4M scores -16 across its junction, so
        // the intronic reference is checked - but reading contiguously (chr1:76-79) also scores -16, which is not higher
        // than the overhang, so the junction is kept.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "20M50N5M60N4M", bases("C".repeat(29)));

        assertEquals("20M50N5M60N4M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testKeptTerminalJunctionLeavesShortInteriorExonUntouched()
    {
        // "20M 50N 4M 60N 40M": the terminal 40M is longer than MIN_OVERHANG_LENGTH so it is trusted unscored and
        // kept, stopping the collapse - the short interior 4M exon flanked by two surviving junctions is never gated.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "20M50N4M60N40M", bases("A".repeat(64)));

        assertEquals("20M50N4M60N40M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testLeadingJunctionCollapses()
    {
        // "4S 5M 100N 20M": leading anchor 5M scores +5 (not > 5) so the 100N collapses; the read-through region
        // (chr1:97-105) is 'C' vs the 'A' read so nothing extends -> "9S 20M" starting at the near exon.
        TestGenome genome = genome().set(CHR_1, 97, "CCCCCCCCC");
        OverhangGate.Result result = gate(genome).collapseJunctions(CHR_1, 1, "4S5M100N20M", bases("A".repeat(29)));

        assertEquals("9S20M", result.cigar());
        assertEquals(106, result.pos());
        assertTrue(result.dropped());
    }

    @Test
    public void testPreCollapsedResolveShape()
    {
        // exp8 read 25535 shape "100M 83N 3M 48S": the 3bp anchor collapses (score < 5), folding into the clip ->
        // "100M 51S".
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1000, "100M83N3M48S", bases("C".repeat(151)));

        assertEquals("100M51S", result.cigar());
        assertTrue(result.dropped());
    }

    @Test
    public void testExtendsThroughMatchingSoftclip()
    {
        // "50M 10S" whose clip continues contiguously in the genome (all 'A' vs all 'A') walks to "60M".
        OverhangGate.Result result = gate(genome()).extendSoftClips(CHR_1, 1, "50M10S", bases("A".repeat(60)));

        assertEquals("60M", result.cigar());
    }

    @Test
    public void testLeavesMismatchedSoftclipUnextended()
    {
        // "50M 10S" whose clip bases ('C') do not match the contiguous genome ('A') is left clipped.
        OverhangGate.Result result = gate(genome())
                .extendSoftClips(CHR_1, 1, "50M10S", bases("A".repeat(50) + "C".repeat(10)));

        assertEquals("50M10S", result.cigar());
    }

    @Test
    public void testDisabledWithoutRefSourceIsNoOp()
    {
        OverhangGate gate = new OverhangGate(null);
        assertFalse(gate.enabled());

        OverhangGate.Result result = gate.collapseJunctions(CHR_1, 1, "20M100N5M5S", bases("C".repeat(30)));
        assertEquals("20M100N5M5S", result.cigar());
        assertFalse(result.dropped());
    }

    // ---- bwa-mem scoring primitives (MATCH +1, MISMATCH -4) ----
    @Test
    public void testScoreSumsMatchesAndMismatchesOverOverlap()
    {
        // 9 matches (+9) then 1 trailing mismatch (-4) over the 10-base overlap sum to +5.
        int result = OverhangGate.score(bases("A".repeat(10)), bases("A".repeat(9) + "C"));
        assertEquals(5, result);
    }

    @Test
    public void testScoreUsesOverlapLengthWhenInputsDiffer()
    {
        // min(10, 4) = 4-base overlap: the 4 leading matches score +4 and the 6 unpaired read bases are not scored.
        int result = OverhangGate.score(bases("A".repeat(10)), bases("AAAA"));
        assertEquals(4, result);
    }

    @Test
    public void testScoreTreatsBasesCaseInsensitively()
    {
        // lowercase read vs uppercase ref: basesEqual counts each as a match (+4).
        int result = OverhangGate.score(bases("aaaa"), bases("AAAA"));
        assertEquals(4, result);
    }

    @Test
    public void testMaxScoringPrefixReturnsFullLengthWhenAllMatch()
    {
        // every base matches, so the whole 10-base window is the highest-scoring prefix.
        int result = OverhangGate.maxScoringPrefix(bases("A".repeat(10)), bases("A".repeat(10)));
        assertEquals(10, result);
    }

    @Test
    public void testMaxScoringPrefixIsZeroWhenNoPrefixScoresPositive()
    {
        // every base mismatches, so the cumulative score never rises above 0 and nothing is extended.
        int result = OverhangGate.maxScoringPrefix(bases("C".repeat(5)), bases("A".repeat(5)));
        assertEquals(0, result);
    }

    @Test
    public void testMaxScoringPrefixStopsAtLongestPositivePrefix()
    {
        // 5 matches (peak +5) then 3 mismatches (dropping to +1, -3, -7): the prefix at the peak is the leading 5.
        int result = OverhangGate.maxScoringPrefix(bases("A".repeat(5) + "C".repeat(3)), bases("A".repeat(8)));
        assertEquals(5, result);
    }

    @Test
    public void testMaxScoringPrefixExtendsMatchesAfterInternalMismatch()
    {
        // 5 matches, 1 mismatch, 4 matches vs all-A: the >= tie climbs back to the +5 peak so the prefix extends to the full 10.
        int result = OverhangGate.maxScoringPrefix(bases("A".repeat(5) + "C" + "A".repeat(4)), bases("A".repeat(10)));
        assertEquals(10, result);
    }

    // ---- genomicScore: bwa-mem score of a lifted placement against the genome ----

    private static int scoreOf(final String cigar, final String sequence)
    {
        return new OverhangGate(genome().asRefGenome()).genomicScore(CHR_1, 1, cigar, bases(sequence));
    }

    @Test
    public void testScoreFullMatchOnePerBase()
    {
        assertEquals(10, scoreOf("10M", "AAAAAAAAAA"));
    }

    @Test
    public void testScoreMismatchCostsFour()
    {
        // one C against all-A ref: 9 matches (+9), 1 mismatch (-4) = 5
        assertEquals(5, scoreOf("10M", "AAAACAAAAA"));
    }

    @Test
    public void testScoreDeletionIsAffineGap()
    {
        // 10 matched bases (+10) and a 2bp deletion (GAP_OPEN -6, GAP_EXTEND -1) = 3
        assertEquals(3, scoreOf("5M2D5M", "AAAAAAAAAA"));
    }

    @Test
    public void testScoreIntronIsFree()
    {
        assertEquals(10, scoreOf("5M100N5M", "AAAAAAAAAA"));
    }

    @Test
    public void testScoreSoftClipIsFree()
    {
        assertEquals(7, scoreOf("3S7M", "AAAAAAAAAA"));
    }

    @Test
    public void testScoreNoRefSourceReturnsSentinel()
    {
        assertEquals(Integer.MIN_VALUE, new OverhangGate(null).genomicScore(CHR_1, 1, "10M", bases("AAAAAAAAAA")));
    }

    // ---- leading side: the mirror of the trailing cases above ----

    @Test
    public void testKeepsLeadingAnchorScoringAboveFive()
    {
        // "4S 6M 100N 20M": the 6bp anchor sits at chr1:1-6 and matches it (score +6 > 5), so the junction is kept.
        TestGenome genome = genome().set(CHR_1, 1, "CCCCCC");
        OverhangGate.Result result = gate(genome).collapseJunctions(CHR_1, 1, "4S6M100N20M", bases("C".repeat(30)));

        assertEquals("4S6M100N20M", result.cigar());
        assertEquals(1, result.pos());
        assertFalse(result.dropped());
    }

    @Test
    public void testCollapsesLeadingShortAnchorOnScore()
    {
        // "4S 8M 100N 20M": an 8bp anchor is short (<= MIN_OVERHANG_LENGTH) so it is scored, mismatches its own exon
        // (all 'A' genome vs 'C' read) and collapses into the clip.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "4S8M100N20M", bases("C".repeat(32)));

        assertEquals("12S20M", result.cigar());
        assertEquals(109, result.pos());
        assertTrue(result.dropped());
    }

    @Test
    public void testTrustsLongLeadingOverhangWithoutScoring()
    {
        // "4S 15M 100N 20M": 15bp exceeds MIN_OVERHANG_LENGTH, so the anchor is trusted unscored and kept even though
        // those bases mismatch the reference.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "4S15M100N20M", bases("C".repeat(39)));

        assertEquals("4S15M100N20M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testLeavesSingleLeadingJunctionNoSoftclipUntouched()
    {
        // "4M 100N 20M": a single junction with no adjacent softclip never triggers the gate.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "4M100N20M", bases("C".repeat(24)));

        assertEquals("4M100N20M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void testLeadingCollapseExtendsIntoTheReadThrough()
    {
        // "4S 5M 100N 20M" with the read-through region matching the read: the collapse extends over those bases as M
        // rather than leaving them clipped, so the clip shrinks and the start moves back.
        OverhangGate.Result result = gate(genome()).collapseJunctions(CHR_1, 1, "4S5M100N20M", bases("A".repeat(29)));

        assertEquals("29M", result.cigar());
        assertEquals(97, result.pos());
        assertTrue(result.dropped());
    }

    @Test
    public void testLeadingCollapseKeepsTrailingElements()
    {
        // "4S 5M 50N 10M 60N 20M": collapsing the first junction must leave the later intron and exon untouched.
        TestGenome genome = genome().set(CHR_1, 47, "CCCCCCCCC");
        OverhangGate.Result result = gate(genome).collapseJunctions(CHR_1, 1, "4S5M50N10M60N20M", bases("A".repeat(39)));

        assertEquals("9S10M60N20M", result.cigar());
        assertEquals(56, result.pos());
        // dropped means every junction went; the 60N survives, so this placement still carries one
        assertFalse(result.dropped());
    }
}
