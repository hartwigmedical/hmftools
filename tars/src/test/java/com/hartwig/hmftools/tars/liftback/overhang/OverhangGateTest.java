package com.hartwig.hmftools.tars.liftback.overhang;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;
import com.hartwig.hmftools.tars.liftback.supplementary.RefSequenceSource;

import org.junit.Test;

// Unit tests for the overhang gate: the absolute anchor-score keep condition (> MIN_OVERHANG_SCORE), the trigger
// (softclip present or multiple junctions), the iterative multi-junction peel, and the standalone softclip reclaim.
public class OverhangGateTest
{
    private static OverhangGate gate(final TestGenome genome)
    {
        return new OverhangGate(genome.asRefSource());
    }

    private static TestGenome genome()
    {
        return new TestGenome().with(CHR_1, 5000, 'A');
    }

    @Test
    public void keepsAnchorScoringAboveFive()
    {
        // trailing "20M 100N 6M 4S": the 6bp anchor matches its far exon (score +6 > 5), so the junction is kept.
        // near exon chr1:1-20, intron 21-120, far exon anchor at chr1:121-126.
        TestGenome genome = genome().set(CHR_1, 121, "CCCCCC");
        OverhangGate.Result result = gate(genome).peel(CHR_1, 1, "20M100N6M4S", bases("C".repeat(30)));

        assertEquals("20M100N6M4S", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void collapsesAnchorScoringExactlyFive()
    {
        // trailing "20M 100N 5M 5S": the 5bp anchor matches its far exon (score +5), which does NOT clear the
        // strict > 5 bar, so the junction collapses; the intron read-through is all 'A' vs 'C' so nothing reclaims.
        TestGenome genome = genome().set(CHR_1, 121, "CCCCC");
        OverhangGate.Result result = gate(genome).peel(CHR_1, 1, "20M100N5M5S", bases("C".repeat(30)));

        assertEquals("20M10S", result.cigar());
        assertTrue(result.dropped());
    }

    @Test
    public void collapsesShortAnchorOnScore()
    {
        // trailing "20M 100N 8M 4S": an 8bp anchor is still short (<= MIN_OVERHANG_LENGTH of 12) so it is scored - it
        // mismatches its far exon (all 'A' vs 'C') so it collapses. Length alone does not save a short anchor.
        OverhangGate.Result result = gate(genome()).peel(CHR_1, 1, "20M100N8M4S", bases("C".repeat(32)));

        assertEquals("20M12S", result.cigar());
        assertTrue(result.dropped());
    }

    @Test
    public void trustsLongOverhangWithoutScoring()
    {
        // trailing "20M 100N 15M 4S": a 15bp anchor is longer than MIN_OVERHANG_LENGTH (12), so it is trusted outright
        // and the junction is kept even though those bases mismatch the far exon (all 'A' vs 'C', which would score < 5).
        OverhangGate.Result result = gate(genome()).peel(CHR_1, 1, "20M100N15M4S", bases("C".repeat(39)));

        assertEquals("20M100N15M4S", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void leavesSingleJunctionNoSoftclipUntouched()
    {
        // "20M 100N 4M": a single junction with no adjacent softclip is Case 3 - left untouched even though the
        // 4bp anchor would fail the score test if it were gated.
        OverhangGate.Result result = gate(genome()).peel(CHR_1, 1, "20M100N4M", bases("C".repeat(24)));

        assertEquals("20M100N4M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void collapsesMultiJunctionWhenReferenceScoresBetter()
    {
        // "20M 50N 5M 60N 4M" (no soft clip, two junctions trigger the gate). The read is all 'A'; the terminal 4M's
        // far exon (chr1:136-139) is set to 'C' so the overhang scores -16, while reading contiguously through the 60N
        // (chr1:76-79) matches (+4). The intronic reference scores higher than the overhang, so the 60N collapses: the
        // 4M reads into the 5M -> "20M 50N 9M". The remaining single junction has no soft clip, so the peel stops.
        TestGenome genome = genome().set(CHR_1, 136, "CCCC");
        OverhangGate.Result result = gate(genome).peel(CHR_1, 1, "20M50N5M60N4M", bases("A".repeat(29)));

        assertEquals("20M50N9M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void keepsMultiJunctionWhenOverhangScoresPositive()
    {
        // "20M 50N 5M 60N 4M", read and genome all 'A': the terminal 4M aligns positively across its junction (+4 > 0),
        // so the junction is kept without even checking the reference alternative.
        OverhangGate.Result result = gate(genome()).peel(CHR_1, 1, "20M50N5M60N4M", bases("A".repeat(29)));

        assertEquals("20M50N5M60N4M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void keepsMultiJunctionWhenReferenceNotBetter()
    {
        // "20M 50N 5M 60N 4M", read all 'C' vs an all-'A' genome: the terminal 4M scores -16 across its junction, so
        // the intronic reference is checked - but reading contiguously (chr1:76-79) also scores -16, which is not higher
        // than the overhang, so the junction is kept.
        OverhangGate.Result result = gate(genome()).peel(CHR_1, 1, "20M50N5M60N4M", bases("C".repeat(29)));

        assertEquals("20M50N5M60N4M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void keptTerminalJunctionLeavesShortInteriorExonUntouched()
    {
        // "20M 50N 4M 60N 40M": the terminal 40M matches (all 'A', score +40) and is kept, so the peel stops
        // immediately - the short interior 4M exon flanked by two surviving junctions is never gated.
        OverhangGate.Result result = gate(genome()).peel(CHR_1, 1, "20M50N4M60N40M", bases("A".repeat(64)));

        assertEquals("20M50N4M60N40M", result.cigar());
        assertFalse(result.dropped());
    }

    @Test
    public void leadingJunctionCollapses()
    {
        // "4S 5M 100N 20M": leading anchor 5M scores +5 (not > 5) so the 100N collapses; the read-through region
        // (chr1:97-105) is 'C' vs the 'A' read so nothing reclaims -> "9S 20M" starting at the near exon.
        TestGenome genome = genome().set(CHR_1, 97, "CCCCCCCCC");
        OverhangGate.Result result = gate(genome).peel(CHR_1, 1, "4S5M100N20M", bases("A".repeat(29)));

        assertEquals("9S20M", result.cigar());
        assertEquals(106, result.pos());
        assertTrue(result.dropped());
    }

    @Test
    public void preCollapsedResolveShape()
    {
        // exp8 read 25535 shape "100M 83N 3M 48S": the 3bp anchor collapses (score < 5), folding into the clip ->
        // "100M 51S", supplementary resolve then merges across the true junction (see SupplementaryResolverTest).
        OverhangGate.Result result = gate(genome()).peel(CHR_1, 1000, "100M83N3M48S", bases("C".repeat(151)));

        assertEquals("100M51S", result.cigar());
        assertTrue(result.dropped());
    }

    @Test
    public void reclaimConvertsMatchingSoftclipPrefix()
    {
        // "50M 10S" whose clip continues contiguously in the genome (all 'A' vs all 'A') walks to "60M".
        OverhangGate.Result result = gate(genome()).reclaimTerminalSoftClip(CHR_1, 1, "50M10S", bases("A".repeat(60)));

        assertEquals("60M", result.cigar());
    }

    @Test
    public void reclaimLeavesMismatchedSoftclip()
    {
        // "50M 10S" whose clip bases ('C') do not match the contiguous genome ('A') is left clipped.
        OverhangGate.Result result = gate(genome())
                .reclaimTerminalSoftClip(CHR_1, 1, "50M10S", bases("A".repeat(50) + "C".repeat(10)));

        assertEquals("50M10S", result.cigar());
    }

    @Test
    public void disabledWithoutRefSourceIsNoOp()
    {
        RefSequenceSource nullSource = null;
        OverhangGate gate = new OverhangGate(nullSource);
        assertFalse(gate.enabled());

        OverhangGate.Result result = gate.peel(CHR_1, 1, "20M100N5M5S", bases("C".repeat(30)));
        assertEquals("20M100N5M5S", result.cigar());
        assertFalse(result.dropped());
    }
}
