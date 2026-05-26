package com.hartwig.hmftools.redux.splice.rescue;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

// Exhaustive coverage for JunctionRescueResolver covering each gate in the design doc edge-case
// table. Reads are 151bp throughout (matches Illumina default for our exp7 data set) unless a
// specific test needs a smaller read for clarity. Tests construct annotated-junction sets inline.
public class JunctionRescueResolverTest
{
    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";
    private static final int READ_LEN = 151;

    private static Set<ChrIntron> annotated(final ChrIntron... introns)
    {
        return new HashSet<>(Arrays.asList(introns));
    }

    private static RescueCandidate candidate(
            final String chrom, final boolean forward, final int readLen, final int primStart,
            final String primCigar, final int primMapq, final RescueSupplementary... supps)
    {
        return new RescueCandidate(chrom, forward, readLen, primStart, primCigar, primMapq,
                supps.length == 0 ? Collections.emptyList() : Arrays.asList(supps));
    }

    private static RescueSupplementary supp(
            final int index, final String chrom, final boolean forward, final int start,
            final String cigar, final int mapq)
    {
        return new RescueSupplementary(index, chrom, forward, start, cigar, mapq);
    }

    private JunctionRescueResolver disabledResolver()
    {
        return new JunctionRescueResolver(Collections.emptySet(), RescueConfig.defaults());
    }

    private JunctionRescueResolver enabledResolver(final Set<ChrIntron> annotated)
    {
        return new JunctionRescueResolver(annotated, RescueConfig.enabledDefaults());
    }

    // ========== happy paths ==========

    @Test
    public void testRightExtendCleanMerge()
    {
        // primary 94M57S + supp 94S57M across an annotated intron; cleanly complementary on read
        // (primary covers read[0..94), supp covers read[94..151)). Tests the core right-extend.
        final int primStart = 31448368;
        final String primCigar = "94M57S";
        final int suppStart = 31448541;
        final String suppCigar = "94S57M";

        // primary ends at primStart + 94 - 1 = 31448461; intron [31448462, 31448540], length 79
        final ChrIntron annotatedIntron = new ChrIntron(CHR1, 31448462, 31448540);
        final RescueSupplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M79N57M", res.MergedCigar);
        assertEquals(primStart, res.MergedStart);
        assertEquals(1, res.DroppedSupplementaryIndices.size());
        assertEquals(Integer.valueOf(0), res.DroppedSupplementaryIndices.get(0));
        assertEquals(1, res.IntroducedIntrons.size());
        assertEquals(annotatedIntron, res.IntroducedIntrons.get(0));
        assertEquals(1, res.ChainDepth);
    }

    @Test
    public void testLeftExtendCleanMerge()
    {
        // mirror of right-extend with complementary cigars: primary 57S94M starts after a junction;
        // supp 57M94S sits upstream covering read[0..57).
        final int suppStart = 31448368;
        final String suppCigar = "57M94S";
        final int primStart = 31448541;
        final String primCigar = "57S94M";

        // supp ref end = 31448368 + 57 - 1 = 31448424; primStart = 31448541
        // intron [31448425, 31448540], length 116
        final ChrIntron intron = new ChrIntron(CHR1, 31448425, 31448540);
        final RescueSupplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("57M116N94M", res.MergedCigar);
        assertEquals(suppStart, res.MergedStart);
    }

    @Test
    public void testChainMergeAcrossThreeExons()
    {
        // primary 50M101S (matches first 50 bases) + supp1 50S60M41S (middle 60 bases) +
        // supp2 110S41M (final 41 bases). Two annotated introns wire it together.
        final int primStart = 1000;
        final String primCigar = "50M101S";
        final int suppMidStart = 2000;        // intron 1: 1050-1999 (len 950)
        final String suppMidCigar = "50S60M41S";
        final int suppLastStart = 3000;       // intron 2: 2060-2999 (len 940)
        final String suppLastCigar = "110S41M";

        final Set<ChrIntron> set = annotated(
                new ChrIntron(CHR1, 1050, 1999),
                new ChrIntron(CHR1, 2060, 2999));

        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255,
                supp(0, CHR1, true, suppMidStart, suppMidCigar, 60),
                supp(1, CHR1, true, suppLastStart, suppLastCigar, 60));

        final RescueResult res = enabledResolver(set).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("50M950N60M940N41M", res.MergedCigar);
        assertEquals(primStart, res.MergedStart);
        assertEquals(2, res.ChainDepth);
        assertEquals(2, res.DroppedSupplementaryIndices.size());
        assertTrue(res.DroppedSupplementaryIndices.contains(0));
        assertTrue(res.DroppedSupplementaryIndices.contains(1));
    }

    @Test
    public void testMergeWhenPrimaryAlreadyHasInternalN()
    {
        // primary already has an internal junction: 50M200N40M61S. The supp picks up the next exon.
        final int primStart = 1000;
        final String primCigar = "50M200N40M61S";
        // primary ref end = 1000 + 50 + 200 + 40 - 1 = 1289
        final int suppStart = 1500;        // intron: 1290-1499 (len 210)
        final String suppCigar = "90S61M";

        final Set<ChrIntron> set = annotated(new ChrIntron(CHR1, 1290, 1499));

        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255,
                supp(0, CHR1, true, suppStart, suppCigar, 60));

        final RescueResult res = enabledResolver(set).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("50M200N40M210N61M", res.MergedCigar);
    }

    // ========== single-gate rejections ==========

    @Test
    public void testRejectWhenDisabled()
    {
        // disabled config returns no merge regardless of inputs
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1175, "90S61M", 60));

        final RescueResult res = disabledResolver().resolve(cand);

        assertFalse(res.Merged);
    }

    @Test
    public void testRejectNoTerminalSoftclip()
    {
        // primary is fully matched (151M); no softclip to extend across
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "151M", 60,
                supp(0, CHR1, true, 2000, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NO_TERMINAL_SOFTCLIP, res.RejectReason);
    }

    @Test
    public void testRejectDifferentChromosome()
    {
        final ChrIntron annotatedIntron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR2, true, 1500, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.DIFFERENT_CHROMOSOME, res.RejectReason);
    }

    @Test
    public void testRejectOppositeStrand()
    {
        final ChrIntron annotatedIntron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, false, 1500, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.OPPOSITE_STRAND, res.RejectReason);
    }

    @Test
    public void testRejectIntronTooShort()
    {
        // intron length 6 — below default MinIntronLength=21. Primary ends at 1093; supp at 1100
        // → intron [1094, 1099].
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1100, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.INTRON_TOO_SHORT, res.RejectReason);
    }

    @Test
    public void testRejectIntronTooLong()
    {
        // intron length 2_000_001 — above default MaxIntronLength=1_000_000
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 2_001_095, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.INTRON_TOO_LONG, res.RejectReason);
    }

    @Test
    public void testRejectShortPrimaryAnchor()
    {
        // primary 2M149S — anchor adjacent to trailing S is only 2 (< MinAnchorOverhang=3)
        final ChrIntron intron = new ChrIntron(CHR1, 1002, 1099);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "2M149S", 60,
                supp(0, CHR1, true, 1100, "2S149M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.SHORT_ANCHOR, res.RejectReason);
    }

    @Test
    public void testRejectShortSuppAnchor()
    {
        // primary anchor 49M is fine; supp leading-M anchor is 1 (< MinAnchorOverhang=3)
        final ChrIntron intron = new ChrIntron(CHR1, 1050, 1099);
        final RescueCandidate cand = candidate(CHR1, true, 50, 1000, "49M1S", 60,
                supp(0, CHR1, true, 1100, "49S1M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.SHORT_ANCHOR, res.RejectReason);
    }

    @Test
    public void testRejectNovelJunctionWhenAnnotatedOnly()
    {
        // Custom config with AnnotatedOnly=true rejects when no L lands on annotated.
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NOVEL_JUNCTION, res.RejectReason);
    }

    @Test
    public void testAcceptNovelJunctionWhenAnnotatedOnlyFalse()
    {
        // disable AnnotatedOnly via custom config; merge should proceed despite empty annotation set.
        final RescueConfig perm = new RescueConfig(true, 21, 1_000_000, 3, 4, false, 0);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), perm).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M106N57M", res.MergedCigar);
    }

    @Test
    public void testRejectHardClipPrimary()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "10H94M57S", 60,
                supp(0, CHR1, true, 1200, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.COMPLEX_CIGAR_SHAPE, res.RejectReason);
    }

    @Test
    public void testRejectHardClipSupp()
    {
        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1199);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1200, "5H90S61M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.COMPLEX_CIGAR_SHAPE, res.RejectReason);
    }

    @Test
    public void testRejectIndelAdjacentToPrimarySoftclip()
    {
        // primary 90M4I57S — has I adjacent to trailing S
        final ChrIntron intron = new ChrIntron(CHR1, 1091, 1199);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "90M4I57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.COMPLEX_CIGAR_SHAPE, res.RejectReason);
    }

    @Test
    public void testRejectReadCoverageOverlap()
    {
        // primary 94M57S covers read[0..94); supp 80S71M covers read[80..151) — overlap 14 bases.
        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "80S71M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.READ_COVERAGE_OVERLAP, res.RejectReason);
    }

    @Test
    public void testRejectReadCoverageGap()
    {
        // primary 94M57S covers read[0..94); supp 110S41M covers read[110..151) — gap 94..110.
        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "110S41M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.READ_COVERAGE_GAP, res.RejectReason);
    }

    @Test
    public void testRejectRefOverlap()
    {
        // primary ends at 1094; supp starts at 1080 — supp upstream of primary's matched end.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1080, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.READ_COVERAGE_OVERLAP, res.RejectReason);
    }

    @Test
    public void testNoSupplementaryAvailable()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60);

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NO_MATCHING_SUPP, res.RejectReason);
    }

    @Test
    public void testRejectShapeMismatchSuppWrongClipSide()
    {
        // primary trailing S expects supp with leading S; this supp has trailing S → invalid pair
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "61M90S", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NO_MATCHING_SUPP, res.RejectReason);
    }

    @Test
    public void testPrimaryBothSidesClippedRightExtendAccepted()
    {
        // Middle-anchored primary in a 3-exon read: both leading and trailing softclip. The right
        // supp's position past primaryRefEnd disambiguates direction → right-extend fires even
        // though the primary has a leading softclip (which the left supp would later pick up).
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                supp(0, CHR1, true, 1500, "95S56M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.Merged);
        // primaryRefEnd = 1090, supp.Start = 1500 → intron length = 1500 - 1 - 1090 = 409.
        assertEquals("5S90M409N56M", res.MergedCigar);
        assertEquals(1001, res.MergedStart);
    }

    @Test
    public void testPrimaryBothSidesClippedLeftExtendAccepted()
    {
        // Same idea on the left. Primary 5S90M56S at chr1:1001 with a left supp ending before
        // primaryStart. Left-extend fires; the trailing-S of primary is the unmerged side.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                supp(0, CHR1, true, 500, "5M146S", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.Merged);
        // suppRefEnd = 504, primaryStart = 1001 → intron length = 1001 - 1 - 504 = 496.
        assertEquals("5M496N90M56S", res.MergedCigar);
        assertEquals(500, res.MergedStart);
    }

    @Test
    public void testPrimaryBothSidesClippedChainMergesBothSupps()
    {
        // Full middle-anchored 3-exon scenario: primary 5S90M56S anchored on the middle exon,
        // right supp at 1500 (95S56M) and left supp at 500 (5M146S). Chain rescue should merge
        // in two iterations: right first, then left. Final cigar spans all 3 exons.
        final RescueSupplementary right = supp(0, CHR1, true, 1500, "95S56M", 60);
        final RescueSupplementary left = supp(1, CHR1, true, 500, "5M146S", 60);
        final RescueCandidate cand = new RescueCandidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                Arrays.asList(right, left));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("5M496N90M409N56M", res.MergedCigar);
        assertEquals(500, res.MergedStart);
        assertEquals(2, res.ChainDepth);
    }

    // ========== ambiguity ==========

    @Test
    public void testAmbiguousSuppChoiceSkipsMerge()
    {
        // two supps with identical MAPQ and identical implied intron length → ambiguous, skip.
        // First supp: start 1500 → intron 1095-1499 (len 405)
        // Second supp: start at a different position producing the same intron length 405:
        //   start 1500 → intron 1095-1499 (same chromosome region) — collide. Use different positions:
        // To make them tie on intron length: use a 2nd supp on a different chromosome? No — must be same chr.
        // Use same chr at start 1500 (intron 405) and start 1500 (same supp data) — but indexes differ.
        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.AMBIGUOUS_SUPP_CHOICE, res.RejectReason);
    }

    @Test
    public void testTwoValidSuppsPicksHighestMapq()
    {
        // supp1 at MAPQ 0, supp2 at MAPQ 60 — resolver should pick supp2.
        final ChrIntron intron1 = new ChrIntron(CHR1, 1095, 1499);
        final ChrIntron intron2 = new ChrIntron(CHR1, 1095, 1799);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 0),
                supp(1, CHR1, true, 1800, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron1, intron2)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals(Integer.valueOf(1), res.DroppedSupplementaryIndices.get(0));
    }

    @Test
    public void testTwoValidSuppsSameMapqPicksShorterIntron()
    {
        // both MAPQ 60; supp1 implies intron 405bp, supp2 implies 805bp — resolver picks supp1.
        final ChrIntron intron1 = new ChrIntron(CHR1, 1095, 1499);
        final ChrIntron intron2 = new ChrIntron(CHR1, 1095, 1899);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1900, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron1, intron2)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals(Integer.valueOf(0), res.DroppedSupplementaryIndices.get(0));
    }

    // ========== zero-MAPQ scenarios (allowed) ==========

    @Test
    public void testMergeWhenPrimaryMapqZero()
    {
        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 0,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.Merged);
    }

    @Test
    public void testMergeWhenSuppMapqZero()
    {
        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 0));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.Merged);
    }

    // ========== chain-depth cap ==========

    @Test
    public void testChainDepthCap()
    {
        // 5 supps would chain into 5-junction read; cap=2 stops the chain after 2 merges.
        final RescueConfig cappedConfig = new RescueConfig(true, 21, 1_000_000, 3, 2, true, 0);

        final int p = 1000;
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, p, "30M121S", 60,
                supp(0, CHR1, true, 2000, "30S30M91S", 60),
                supp(1, CHR1, true, 3000, "60S30M61S", 60),
                supp(2, CHR1, true, 4000, "90S30M31S", 60),
                supp(3, CHR1, true, 5000, "120S31M", 60));

        final Set<ChrIntron> set = annotated(
                new ChrIntron(CHR1, 1030, 1999),
                new ChrIntron(CHR1, 2030, 2999),
                new ChrIntron(CHR1, 3030, 3999),
                new ChrIntron(CHR1, 4030, 4999));

        final RescueResult res = new JunctionRescueResolver(set, cappedConfig).resolve(cand);

        assertTrue(res.Merged);
        assertEquals(2, res.ChainDepth);
        assertEquals(2, res.DroppedSupplementaryIndices.size());
    }

    // ========== other read lengths ==========

    @Test
    public void testShortReadLength()
    {
        // 50bp read with 30M20S + 30S20M across a 100bp intron
        final ChrIntron intron = new ChrIntron(CHR1, 1031, 1130);
        final RescueCandidate cand = candidate(CHR1, true, 50, 1001, "30M20S", 60,
                supp(0, CHR1, true, 1131, "30S20M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("30M100N20M", res.MergedCigar);
    }

    // ========== statistics ==========

    @Test
    public void testStatisticsCounters()
    {
        // strict config (AnnotatedOnly=true) so the novel-junction reject is observable.
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5);
        final JunctionRescueResolver resolver = new JunctionRescueResolver(
                annotated(new ChrIntron(CHR1, 1095, 1499)), strict);

        // success
        resolver.resolve(candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60)));
        // novel-junction reject
        resolver.resolve(candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1800, "94S57M", 60)));
        // diff-chromosome reject
        resolver.resolve(candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR2, true, 1500, "94S57M", 60)));

        final RescueStatistics stats = resolver.statistics();
        assertEquals(3, stats.candidatesEvaluated());
        assertEquals(1, stats.mergedTotal());
        assertEquals(1, stats.mergedAtChainDepth(1));
        assertEquals(1, stats.rejectCount(RescueRejectReason.NOVEL_JUNCTION));
        assertEquals(1, stats.rejectCount(RescueRejectReason.DIFFERENT_CHROMOSOME));
    }

    @Test
    public void testNoMergeStillCountsCandidateAndReject()
    {
        final JunctionRescueResolver resolver = new JunctionRescueResolver(
                Collections.emptySet(), RescueConfig.enabledDefaults());

        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, false, 1500, "94S57M", 60));

        final RescueResult res = resolver.resolve(cand);

        assertFalse(res.Merged);
        assertEquals(1, resolver.statistics().candidatesEvaluated());
        assertEquals(0, resolver.statistics().mergedTotal());
        assertEquals(1, resolver.statistics().rejectCount(RescueRejectReason.OPPOSITE_STRAND));
    }

    // ========== exp7 example reads ==========

    @Test
    public void testExp7Case2Chr1_31448368()
    {
        // Inspired by exp7 case 2 (chr1:31448368). The real BWA output had a 4bp overlap between
        // primary and supp matched regions; v1 doesn't handle overlap (see testExp7Case3 for the
        // explicit overlap-rejection case). This test uses the clean complementary cigars BWA
        // would have emitted with a different seed, demonstrating the right-extend merge produces
        // the spliced cigar STAR ends up with.
        final ChrIntron annotatedIntron = new ChrIntron(CHR1, 31448462, 31448539);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 31448368, "94M57S", 60,
                supp(0, CHR1, true, 31448540, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M78N57M", res.MergedCigar);
        assertEquals(31448368, res.MergedStart);
    }

    @Test
    public void testExp7Case3Chr5_34937631()
    {
        // Real exp7 case: chr5:34937631 read with primary 59S92M at chr5:34938856 + supp 61M90S
        // at chr5:34937631. Primary leading-S=59, supp matched=61, overlap=2. STAR's spliced
        // primary at chr5:34937631 is 61M1166N90M — intron lives at (34937692, 34938857). The
        // snap finds L=61 (trust-supp end of overlap) lands on the annotated junction, so the
        // merge produces STAR's exact CIGAR.
        final String chr5 = "chr5";
        final ChrIntron annotatedIntron = new ChrIntron(chr5, 34937692, 34938857);
        final RescueCandidate cand = candidate(chr5, true, READ_LEN, 34938856, "59S92M", 60,
                supp(0, chr5, true, 34937631, "61M90S", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("61M1166N90M", res.MergedCigar);
        assertEquals(34937631, res.MergedStart);
    }

    // ========== overlap-tolerance behavior ==========

    @Test
    public void testOverlapWithinToleranceSnapsToAnnotatedJunction()
    {
        // primary 94M57S at 1001 + supp 90S61M at 1500. overlap = 4.
        // intronLength = (1500-1-1094) + 4 = 409.
        // L=94 (trust primary): primaryLoss=0, suppLoss=4 → intron (1095, 1503).
        final ChrIntron annotated = new ChrIntron(CHR1, 1095, 1503);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated(annotated)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M409N57M", res.MergedCigar);
        assertEquals(1001, res.MergedStart);
        assertEquals(annotated, res.IntroducedIntrons.get(0));
    }

    @Test
    public void testOverlapWithinToleranceSnapsToAnnotatedAtTrustSuppEnd()
    {
        // primary 94M57S at 1000 + supp 92S59M at 1498. overlap = 2.
        // intronLength = (1498-1-1093) + 2 = 406.
        // L=94 (trust primary): intron (1095, 1500).
        // L=92 (trust supp): primaryLoss=2, suppLoss=0 → intron (1092, 1497).
        final ChrIntron annotated = new ChrIntron(CHR1, 1092, 1497);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = enabledResolver(annotated(annotated)).resolve(cand);

        assertTrue(res.Merged);
        // primaryLoss=2 → primary 94M → 92M; suppLoss=0 → supp 59M unchanged.
        assertEquals("92M406N59M", res.MergedCigar);
    }

    @Test
    public void testOverlapWithinToleranceMaxMinAnchorWins()
    {
        // primary 94M57S at 1001 + supp 92S59M at 1498. overlap = 2.
        // intronLength = (1498-1-1094) + 2 = 405.
        // L=94 (trust primary): primary anchor 94, supp anchor 57-2=55. min = 55.
        // L=92 (trust supp): primary anchor 94-2=92, supp anchor 59. min = 59.
        // Both annotated → max-min-anchor rule picks trust-supp (L=92, min anchor 59).
        // Equivalent to STAR: the split point with the strongest support on both sides.
        final ChrIntron trustPrimary = new ChrIntron(CHR1, 1095, 1499);
        final ChrIntron trustSupp = new ChrIntron(CHR1, 1093, 1497);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = enabledResolver(annotated(trustPrimary, trustSupp)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("92M405N59M", res.MergedCigar);
        assertEquals(trustSupp, res.IntroducedIntrons.get(0));
    }

    @Test
    public void testOverlapExceedsToleranceRejected()
    {
        // 6-base overlap exceeds default tolerance (5).
        final ChrIntron annotated = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "88S63M", 60));

        final RescueResult res = enabledResolver(annotated(annotated)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.READ_COVERAGE_OVERLAP, res.RejectReason);
    }

    @Test
    public void testNoAnnotatedSnapWithAnnotatedOnlyTrueRejects()
    {
        // Strict config (AnnotatedOnly=true) rejects when no L lands on annotated.
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NOVEL_JUNCTION, res.RejectReason);
    }

    @Test
    public void testNoAnnotatedSnapWithDefaultsFallsBackToTrustPrimary()
    {
        // Default enabled config has AnnotatedOnly=false — merge succeeds via trust-primary
        // fallback when no L lands on annotated.
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.Merged);
        // trust-primary fallback: primaryLoss=0, suppLoss=2 → 94M + N(407) + 57M
        assertEquals("94M407N57M", res.MergedCigar);
    }

    @Test
    public void testNoAnnotatedSnapWithAnnotatedOnlyFalseFallsBackToTrustPrimary()
    {
        // AnnotatedOnly=false → if no L is annotated, fall back to trust-primary (L=primaryMatched).
        // primary 94M57S at 1001 + supp 92S59M at 1500. overlap = 2.
        // intronLength = (1500-1-1094) + 2 = 407.
        final RescueConfig perm = new RescueConfig(true, 21, 1_000_000, 3, 4, false, 5);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), perm).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M407N57M", res.MergedCigar);
    }

    // ========== misc ==========

    @Test
    public void testMergedIntronListReflectsChainOrder()
    {
        final int primStart = 1000;
        final Set<ChrIntron> set = annotated(
                new ChrIntron(CHR1, 1050, 1999),
                new ChrIntron(CHR1, 2060, 2999));

        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, "50M101S", 60,
                supp(0, CHR1, true, 2000, "50S60M41S", 60),
                supp(1, CHR1, true, 3000, "110S41M", 60));

        final RescueResult res = enabledResolver(set).resolve(cand);

        assertTrue(res.Merged);
        assertEquals(2, res.IntroducedIntrons.size());
        assertEquals(new ChrIntron(CHR1, 1050, 1999), res.IntroducedIntrons.get(0));
        assertEquals(new ChrIntron(CHR1, 2060, 2999), res.IntroducedIntrons.get(1));
    }

    // ========== Layer 1: ref-verify path ==========

    private static RefSequenceSource refSource(final java.util.Map<String, byte[]> chromBases)
    {
        return (chrom, posStart, posEnd) ->
        {
            final byte[] bases = chromBases.get(chrom);
            if(bases == null || posStart < 1 || posEnd > bases.length)
                return null;
            final byte[] out = new byte[posEnd - posStart + 1];
            System.arraycopy(bases, posStart - 1, out, 0, out.length);
            return out;
        };
    }

    private static byte[] bytes(final String s)
    {
        return s.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
    }

    private static JunctionRescueResolver refVerifyResolver(
            final java.util.Set<ChrIntron> annotated,
            final java.util.Map<String, byte[]> chromBases)
    {
        return new JunctionRescueResolver(
                new AnnotatedJunctionIndex(annotated),
                refSource(chromBases),
                RescueConfig.enabledDefaults());
    }

    // ========== splice-motif scan (cryptic-junction recovery) ==========

    // builds a ref backed by uniform 'N' bases, then writes the canonical GT...AG motif at the
    // intron's flank positions so the motif scan tier fires. Returns the chromosome bases map
    // ready for refSource().
    private static java.util.Map<String, byte[]> refWithCanonicalIntron(
            final int chromLen, final int intronStart, final int intronEnd)
    {
        final byte[] bases = new byte[chromLen];
        java.util.Arrays.fill(bases, (byte) 'N');
        bases[intronStart - 1] = 'G'; bases[intronStart] = 'T';        // donor at intronStart..+1 (1-based)
        bases[intronEnd - 2] = 'A'; bases[intronEnd - 1] = 'G';         // acceptor at intronEnd-1..end
        return java.util.Collections.singletonMap(CHR1, bases);
    }

    private static JunctionRescueResolver motifResolver(
            final java.util.Set<ChrIntron> annotated, final java.util.Map<String, byte[]> chromBases)
    {
        return new JunctionRescueResolver(
                new AnnotatedJunctionIndex(annotated), refSource(chromBases),
                RescueConfig.enabledDefaults());
    }

    @Test
    public void testMotifScanPicksCanonicalGTagWhenUnannotated()
    {
        // primary 94M57S at 1001, supp 94S57M at 1500. overlap=0. intron (1095, 1499) length 405.
        // Ref placed so 1095..1096="GT" and 1498..1499="AG" — canonical motif. Merge should fire
        // via the motif scan path even though the annotated index is empty.
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = motifResolver(
                java.util.Collections.emptySet(), refWithCanonicalIntron(2000, 1095, 1499)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M405N57M", res.MergedCigar);
        assertEquals(new ChrIntron(CHR1, 1095, 1499), res.IntroducedIntrons.get(0));
    }

    @Test
    public void testMotifScanPicksAnnotatedOverMotifWhenBothPresent()
    {
        // overlap=2 with two viable L's: L=94 (trust primary) and L=92 (trust supp). Place an
        // ANNOTATED intron at L=94 and a CANONICAL MOTIF at L=92. Annotation must win regardless
        // of min-anchor.
        final ChrIntron annotated = new ChrIntron(CHR1, 1095, 1499);
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        // motif candidate at (1093, 1497): donor 1093..1094 = "GT", acceptor 1496..1497 = "AG"
        bases[1092] = 'G'; bases[1093] = 'T';
        bases[1495] = 'A'; bases[1496] = 'G';

        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = motifResolver(
                annotated(annotated), java.util.Collections.singletonMap(CHR1, bases)).resolve(cand);

        assertTrue(res.Merged);
        // annotated at L=94 wins over motif at L=92
        assertEquals("94M405N57M", res.MergedCigar);
        assertEquals(annotated, res.IntroducedIntrons.get(0));
    }

    @Test
    public void testMotifScanPrefersCanonicalOverSemiCanonical()
    {
        // overlap=2 with two viable L's. Place a CANONICAL (GT-AG) motif at L=92 and a
        // SEMI-CANONICAL (GC-AG) motif at L=94. Canonical (tier 2) must beat semi-canonical (tier 1)
        // regardless of min-anchor.
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        // semi at (1095, 1499): donor 1095..1096 = "GC", acceptor 1498..1499 = "AG"
        bases[1094] = 'G'; bases[1095] = 'C';
        bases[1497] = 'A'; bases[1498] = 'G';
        // canonical at (1093, 1497): donor 1093..1094 = "GT", acceptor 1496..1497 = "AG"
        bases[1092] = 'G'; bases[1093] = 'T';
        bases[1495] = 'A'; bases[1496] = 'G';

        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = motifResolver(
                java.util.Collections.emptySet(), java.util.Collections.singletonMap(CHR1, bases)).resolve(cand);

        assertTrue(res.Merged);
        // canonical at L=92 wins over semi at L=94
        assertEquals("92M405N59M", res.MergedCigar);
        assertEquals(new ChrIntron(CHR1, 1093, 1497), res.IntroducedIntrons.get(0));
    }

    @Test
    public void testMotifScanAcceptsReverseStrandCanonical()
    {
        // Reverse-strand transcript: genomic forward shows "CT" at donor and "AC" at acceptor
        // (RC of AG and GT). overlap=0, intron (1095, 1499). Merge via motif scan.
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        bases[1094] = 'C'; bases[1095] = 'T';    // intronStart=1095..1096 = "CT"
        bases[1497] = 'A'; bases[1498] = 'C';    // intronEnd=1498..1499 = "AC"

        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = motifResolver(
                java.util.Collections.emptySet(), java.util.Collections.singletonMap(CHR1, bases)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M405N57M", res.MergedCigar);
    }

    @Test
    public void testMotifScanIgnoredWhenNoMotifAndNoAnnotated()
    {
        // overlap=0, ref bases all 'N', no annotation, defaults (AnnotatedOnly=false). Snap loop
        // finds no tier-1+ candidate, falls through to trust-primary fallback. Behaviour unchanged.
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = motifResolver(
                java.util.Collections.emptySet(), java.util.Collections.singletonMap(CHR1, bases)).resolve(cand);

        assertTrue(res.Merged);
        // trust-primary fallback: L=primaryMatched=94, intron (1095, 1499), length 405
        assertEquals("94M405N57M", res.MergedCigar);
    }

    @Test
    public void testMotifScanLeftExtendCanonical()
    {
        // Left-extend with no overlap: primary 57S94M at 1500, supp 57M94S at 1001.
        // suppRefEnd = 1057, primaryStart = 1500 → intron (1058, 1499), length 442.
        // Ref: 1058..1059 = "GT", 1498..1499 = "AG" — canonical motif.
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        bases[1057] = 'G'; bases[1058] = 'T';
        bases[1497] = 'A'; bases[1498] = 'G';

        final RescueCandidate cand = candidate(CHR1, true, 151, 1500, "57S94M", 60,
                supp(0, CHR1, true, 1001, "57M94S", 60));

        final RescueResult res = motifResolver(
                java.util.Collections.emptySet(), java.util.Collections.singletonMap(CHR1, bases)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("57M442N94M", res.MergedCigar);
    }

    @Test
    public void testRefVerifyRightExtendSuccess()
    {
        // Primary 30M5S at chr1:101. Trailing 5S read bases "CCCCC".
        // Annotated intron (chr1, 131, 200). Candidate downstream exon at chr1:201.
        // We supply ref bases such that chr1:201..205 = "CCCCC" exactly.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        // Place "CCCCC" at position 201..205 (0-based 200..204)
        for(int i = 200; i < 205; ++i) chr1Ref[i] = 'C';

        final byte[] readBases = new byte[35];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';     // primary anchor
        for(int i = 30; i < 35; ++i) readBases[i] = 'C';    // trailing softclip

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 35, 101, "30M5S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("30M70N5M", res.MergedCigar);
        assertEquals(101, res.MergedStart);
    }

    @Test
    public void testRefVerifyRightExtendRejectsOnMismatch()
    {
        // Same setup but ref bases at downstream exon are "GGGGG" — read's "CCCCC" can't match.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'G';

        final byte[] readBases = new byte[35];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';
        for(int i = 30; i < 35; ++i) readBases[i] = 'C';

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 35, 101, "30M5S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH, res.RejectReason);
    }

    @Test
    public void testRefVerifyAllows10PercentMismatch()
    {
        // 20-base softclip: tolerance = floor(20/10) = 2 mismatches.
        // Read "CCCCCCCCCCCCCCCCCCCC", ref "CCCCCCCCCCCCCCCCCCGC" — 1 mismatch, accept.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        // ref at 201..220 = "CCCCCCCCCCCCCCCCCCGC" (positions 200..219 0-based)
        for(int i = 200; i < 218; ++i) chr1Ref[i] = 'C';
        chr1Ref[218] = 'G';     // mismatch
        chr1Ref[219] = 'C';

        final byte[] readBases = new byte[40];
        for(int i = 0; i < 20; ++i) readBases[i] = 'A';
        for(int i = 20; i < 40; ++i) readBases[i] = 'C';

        final ChrIntron intron = new ChrIntron(CHR1, 121, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 40, 101, "20M20S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("20M80N20M", res.MergedCigar);
    }

    @Test
    public void testRefVerifyLeftExtendSuccess()
    {
        // Primary 5S30M at chr1:201. Leading 5S read bases "TTTTT".
        // Annotated intron (chr1, 131, 200), upstream exon ends at chr1:130.
        // Place "TTTTT" at chr1:126..130 (0-based 125..129).
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        for(int i = 125; i < 130; ++i) chr1Ref[i] = 'T';

        final byte[] readBases = new byte[35];
        for(int i = 0; i < 5; ++i) readBases[i] = 'T';
        for(int i = 5; i < 35; ++i) readBases[i] = 'A';

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 35, 201, "5S30M", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("5M70N30M", res.MergedCigar);
        assertEquals(126, res.MergedStart);
    }

    @Test
    public void testRefVerifyRejectsWhenNoCandidateExon()
    {
        // Primary 30M5S at chr1:101 — annotation has no intron starting at 131, so no candidate.
        final byte[] readBases = new byte[35];
        final byte[] chr1Ref = new byte[300];
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 35, 101, "30M5S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                Collections.emptySet(),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON, res.RejectReason);
    }

    @Test
    public void testRefVerifyAmbiguousRejected()
    {
        // Two annotated introns share the same start; their downstream exons both match the
        // softclipped read bases equally well → ambiguous → skip.
        final byte[] chr1Ref = new byte[600];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        // Ref bases at 201..205 and at 501..505 both equal "CCCCC"
        for(int i = 200; i < 205; ++i) chr1Ref[i] = 'C';
        for(int i = 500; i < 505; ++i) chr1Ref[i] = 'C';

        final byte[] readBases = new byte[35];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';
        for(int i = 30; i < 35; ++i) readBases[i] = 'C';

        final ChrIntron intron1 = new ChrIntron(CHR1, 131, 200);
        final ChrIntron intron2 = new ChrIntron(CHR1, 131, 500);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 35, 101, "30M5S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron1, intron2)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.REF_VERIFY_AMBIGUOUS, res.RejectReason);
    }

    @Test
    public void testRefVerifyWithoutRefSourceSkipsSilently()
    {
        // No RefSequenceSource configured → ref-verify can't run, return NO_MATCHING_SUPP.
        final RescueCandidate cand = candidate(CHR1, true, 35, 101, "30M5S", 60);

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NO_MATCHING_SUPP, res.RejectReason);
    }

    // ========== Layer 2: mate-informed snap hint ==========

    @Test
    public void testMateHintUsedAsFallbackWhenAnnotationMisses()
    {
        // overlap=2, no annotated junction within window. Without mate hint → trust-primary
        // fallback. With mate hint pointing to the L=92 (trust-supp) split → use the hint.
        final RescueCandidate withoutHint = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));
        final RescueResult resNoHint = enabledResolver(annotated()).resolve(withoutHint);
        assertTrue(resNoHint.Merged);
        // trust-primary fallback: intron starts at 1095
        assertEquals(1095, resNoHint.IntroducedIntrons.get(0).IntronStart);

        // With hint pointing at the trust-supp intron (1093, 1497)
        final ChrIntron hint = new ChrIntron(CHR1, 1093, 1497);
        final RescueCandidate withHint = new RescueCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));
        final RescueResult resHinted = enabledResolver(annotated()).resolve(withHint);
        assertTrue(resHinted.Merged);
        assertEquals(1093, resHinted.IntroducedIntrons.get(0).IntronStart);
        // trust-supp snap: primary anchor 94M → 92M, supp 59M unchanged
        assertEquals("92M405N59M", resHinted.MergedCigar);
    }

    @Test
    public void testMateHintIgnoredWhenAnnotatedMatchExists()
    {
        // Annotated junction wins over mate hint when both are within the snap window.
        final ChrIntron annotatedJunc = new ChrIntron(CHR1, 1095, 1499);
        final ChrIntron hint = new ChrIntron(CHR1, 1093, 1497);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));

        final RescueResult res = enabledResolver(annotated(annotatedJunc)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals(1095, res.IntroducedIntrons.get(0).IntronStart);
    }

    @Test
    public void testMateHintOutsideOverlapWindowIgnored()
    {
        // Hint intron is way outside the read's overlap window → resolver falls back to trust-primary.
        final ChrIntron hint = new ChrIntron(CHR1, 50000, 50100);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.Merged);
        assertEquals(1095, res.IntroducedIntrons.get(0).IntronStart);
    }

    @Test
    public void testRejectReasonIsNullOnSuccess()
    {
        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.Merged);
        assertNull(res.RejectReason);
    }
}
