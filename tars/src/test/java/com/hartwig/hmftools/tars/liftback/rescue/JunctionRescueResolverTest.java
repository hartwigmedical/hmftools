package com.hartwig.hmftools.tars.liftback.rescue;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import com.hartwig.hmftools.tars.common.SpliceCommon;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures;

import org.junit.Test;

// Tests for JunctionRescueResolver. Reads are 151bp (Illumina default) unless noted.
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

    @Test
    public void testRightExtendCleanMerge()
    {
        // Clean complementary cigars across an annotated intron; tests core right-extend.
        final int primStart = 31448368;
        final String primCigar = "94M57S";
        final int suppStart = 31448541;
        final String suppCigar = "94S57M";

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
    public void testFoldTrailingFabricatedMicroJunction()
    {
        // 100M83N3M48S: 3bp anchor split off by a spurious 83N against the trailing clip -> 100M51S.
        assertEquals("100M51S", CigarShape.format(JunctionRescueResolver.foldFabricatedTerminalMicroJunctions(
                CigarShape.parse("100M83N3M48S"), SpliceCommon.MIN_JUNCTION_ANCHOR)));
    }

    @Test
    public void testFoldLeadingFabricatedMicroJunction()
    {
        assertEquals("51S100M", CigarShape.format(JunctionRescueResolver.foldFabricatedTerminalMicroJunctions(
                CigarShape.parse("48S3M83N100M"), SpliceCommon.MIN_JUNCTION_ANCHOR)));
    }

    @Test
    public void testNoFoldWhenAnchorAboveThreshold()
    {
        // 8bp anchor meets MIN_JUNCTION_ANCHOR (8) — a trusted junction, left intact.
        assertEquals("100M83N8M48S", CigarShape.format(JunctionRescueResolver.foldFabricatedTerminalMicroJunctions(
                CigarShape.parse("100M83N8M48S"), SpliceCommon.MIN_JUNCTION_ANCHOR)));
    }

    @Test
    public void testNoFoldWithoutTerminalSoftClip()
    {
        // no terminal softclip — nothing to fold into.
        assertEquals("100M83N3M", CigarShape.format(JunctionRescueResolver.foldFabricatedTerminalMicroJunctions(
                CigarShape.parse("100M83N3M"), SpliceCommon.MIN_JUNCTION_ANCHOR)));
    }

    @Test
    public void testRescueAcrossFabricatedTerminalMicroJunction()
    {
        // exp8 read 25535: tx-contig over-run produces 100M83N3M48S; the 3bp fabricated micro-anchor
        // defeats supp merge. Folding to 100M51S lets the merge find the true 156N junction.
        final int primStart = 1051270;
        final String primCigar = "100M83N3M48S";
        final int suppStart = 1051525;
        final String suppCigar = "99S52M";

        final ChrIntron annotatedIntron = new ChrIntron(CHR1, 1051370, 1051525);
        final RescueSupplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("100M156N51M", res.MergedCigar);
        assertEquals(primStart, res.MergedStart);
        assertEquals(annotatedIntron, res.IntroducedIntrons.get(0));
    }

    @Test
    public void testLeftExtendCleanMerge()
    {
        // Mirror of right-extend: primary starts after junction, supp covers the upstream exon.
        final int suppStart = 31448368;
        final String suppCigar = "57M94S";
        final int primStart = 31448541;
        final String primCigar = "57S94M";

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
        // 3-exon read: primary on first exon, two supps on middle and last exons.
        final int primStart = 1000;
        final String primCigar = "50M101S";
        final int suppMidStart = 2000;
        final String suppMidCigar = "50S60M41S";
        final int suppLastStart = 3000;
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
        // Primary already has an internal junction; supp picks up the next exon.
        final int primStart = 1000;
        final String primCigar = "50M200N40M61S";
        final int suppStart = 1500;
        final String suppCigar = "90S61M";

        final Set<ChrIntron> set = annotated(new ChrIntron(CHR1, 1290, 1499));

        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255,
                supp(0, CHR1, true, suppStart, suppCigar, 60));

        final RescueResult res = enabledResolver(set).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("50M200N40M210N61M", res.MergedCigar);
    }

    @Test
    public void testRejectWhenDisabled()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1175, "90S61M", 60));

        final RescueResult res = disabledResolver().resolve(cand);

        assertFalse(res.Merged);
    }

    @Test
    public void testRejectNoTerminalSoftclip()
    {
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
        // Intron length 6 — below default MinIntronLength=21.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1100, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.INTRON_TOO_SHORT, res.RejectReason);
    }

    @Test
    public void testRejectIntronTooLong()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 2_001_095, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.INTRON_TOO_LONG, res.RejectReason);
    }

    @Test
    public void testRejectShortPrimaryAnchor()
    {
        // Primary anchor 2bp — below MinAnchorOverhang=3.
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
        // Supp anchor 1bp — below MinAnchorOverhang=3.
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
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5, 0, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NOVEL_JUNCTION, res.RejectReason);
    }

    @Test
    public void testAcceptNovelJunctionWhenAnnotatedOnlyFalse()
    {
        final RescueConfig perm = new RescueConfig(true, 21, 1_000_000, 3, 4, false, 0, 0, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
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
        // overlap 14 bases — exceeds tolerance.
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
        // Supp starts upstream of primary's matched end — ref overlap.
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
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "61M90S", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NO_MATCHING_SUPP, res.RejectReason);
    }

    @Test
    public void testPrimaryBothSidesClippedRightExtendAccepted()
    {
        // Middle-anchored primary: supp past primaryRefEnd disambiguates direction → right-extend fires.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                supp(0, CHR1, true, 1500, "95S56M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("5S90M409N56M", res.MergedCigar);
        assertEquals(1001, res.MergedStart);
    }

    @Test
    public void testPrimaryBothSidesClippedLeftExtendAccepted()
    {
        // Mirror on the left: supp ends before primaryStart → left-extend fires.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                supp(0, CHR1, true, 500, "5M146S", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("5M496N90M56S", res.MergedCigar);
        assertEquals(500, res.MergedStart);
    }

    @Test
    public void testPrimaryBothSidesClippedChainMergesBothSupps()
    {
        // Full middle-anchored 3-exon scenario: chain merges right first then left.
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

    @Test
    public void testMultipleSuppsInReachSkipsMerge()
    {
        // Two supps both pass every gate — refuse to guess the splice destination.
        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.MULTIPLE_SUPPS_IN_REACH, res.RejectReason);
    }

    @Test
    public void testPrimaryMapqZeroStillMerges()
    {
        // MAPQ-0 primary is the common tx-contig duplicate artifact — not a disqualifier.
        final ChrIntron intron = new ChrIntron(CHR1, 31448462, 31448540);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 31448368, "94M57S", 0,
                supp(0, CHR1, true, 31448541, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.Merged);
    }

    @Test
    public void testTwoValidSuppsSameBoundaryRefusesToGuess()
    {
        // Two supps at different intron lengths both pass — single-supp-within-reach policy refuses to guess.
        final ChrIntron intron1 = new ChrIntron(CHR1, 1095, 1499);
        final ChrIntron intron2 = new ChrIntron(CHR1, 1095, 1799);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1800, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron1, intron2)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.MULTIPLE_SUPPS_IN_REACH, res.RejectReason);
    }

    @Test
    public void testPrimaryMapqFloorKnobStillRejects()
    {
        // Default floor is 0 (MAPQ-0 primaries merge), but a floor of 1 rejects them.
        final RescueConfig floorConfig = new RescueConfig(
                true, RescueConfig.DEFAULT_MIN_INTRON_LENGTH, RescueConfig.DEFAULT_MAX_INTRON_LENGTH,
                RescueConfig.DEFAULT_MIN_ANCHOR_OVERHANG, RescueConfig.DEFAULT_MAX_CHAIN_DEPTH, false,
                RescueConfig.DEFAULT_SOFTCLIP_TOLERANCE, RescueConfig.DEFAULT_MAX_BOUNDARY_SHIFT, 1,
                RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);

        final ChrIntron intron = new ChrIntron(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 0,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = new JunctionRescueResolver(annotated(intron), floorConfig).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.LOW_PRIMARY_MAPQ, res.RejectReason);
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

    @Test
    public void testChainDepthCap()
    {
        // cap=2 stops the chain after 2 merges even when more supps are available.
        final RescueConfig cappedConfig = new RescueConfig(true, 21, 1_000_000, 3, 2, true, 0, 0, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);

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

    @Test
    public void testShortReadLength()
    {
        final ChrIntron intron = new ChrIntron(CHR1, 1031, 1130);
        final RescueCandidate cand = candidate(CHR1, true, 50, 1001, "30M20S", 60,
                supp(0, CHR1, true, 1131, "30S20M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("30M100N20M", res.MergedCigar);
    }

    @Test
    public void testStatisticsCounters()
    {
        // AnnotatedOnly=true so the novel-junction reject is observable.
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5, 0, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final JunctionRescueResolver resolver = new JunctionRescueResolver(
                annotated(new ChrIntron(CHR1, 1095, 1499)), strict);

        resolver.resolve(candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60)));
        resolver.resolve(candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1800, "94S57M", 60)));
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

    @Test
    public void testExp7Case2Chr1_31448368()
    {
        // exp7 case 2 (chr1:31448368): clean complementary cigars produce the same CIGAR as STAR.
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
        // exp7 chr5:34937631: overlap=2, snap picks L=61 (trust-supp) on the annotated junction,
        // producing STAR's exact 61M1166N90M.
        final String chr5 = "chr5";
        final ChrIntron annotatedIntron = new ChrIntron(chr5, 34937692, 34938857);
        final RescueCandidate cand = candidate(chr5, true, READ_LEN, 34938856, "59S92M", 60,
                supp(0, chr5, true, 34937631, "61M90S", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("61M1166N90M", res.MergedCigar);
        assertEquals(34937631, res.MergedStart);
    }

    @Test
    public void testOverlapWithinToleranceSnapsToAnnotatedJunction()
    {
        // overlap=4; trust-primary L=94 lands on annotated (1095, 1503).
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
        // overlap=2; trust-supp L=92 lands on annotated (1092, 1497).
        final ChrIntron annotated = new ChrIntron(CHR1, 1092, 1497);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = enabledResolver(annotated(annotated)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("92M406N59M", res.MergedCigar);
    }

    @Test
    public void testOverlapWithinToleranceMaxMinAnchorWins()
    {
        // overlap=2, both L's annotated. max-min-anchor: trust-supp (min=59) beats trust-primary (min=55).
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
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5, 0, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NOVEL_JUNCTION, res.RejectReason);
    }

    @Test
    public void testNoAnnotatedSnapWithDefaultsFallsBackToTrustPrimary()
    {
        // AnnotatedOnly=false — merge succeeds via trust-primary fallback when no L is annotated.
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M407N57M", res.MergedCigar);
    }

    @Test
    public void testNoAnnotatedSnapWithAnnotatedOnlyFalseFallsBackToTrustPrimary()
    {
        final RescueConfig perm = new RescueConfig(true, 21, 1_000_000, 3, 4, false, 5, 0, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), perm).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M407N57M", res.MergedCigar);
    }

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

    // single-chromosome in-memory ref backed by the shared fixture (MockRefGenome, 1-based inclusive).
    private static RefSequenceSource refSource(final java.util.Map<String, byte[]> chromBases)
    {
        final java.util.Map.Entry<String, byte[]> entry = chromBases.entrySet().iterator().next();
        return TarsTestFixtures.refSource(entry.getKey(), new String(entry.getValue(), java.nio.charset.StandardCharsets.US_ASCII));
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

    // Builds a ref of 'N' bases with GT...AG canonical motif at the intron flanks.
    private static java.util.Map<String, byte[]> refWithCanonicalIntron(
            final int chromLen, final int intronStart, final int intronEnd)
    {
        final byte[] bases = new byte[chromLen];
        java.util.Arrays.fill(bases, (byte) 'N');
        bases[intronStart - 1] = 'G'; bases[intronStart] = 'T';
        bases[intronEnd - 2] = 'A'; bases[intronEnd - 1] = 'G';
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
        // Canonical GT-AG motif placed at intron (1095, 1499) with no annotation — merge via motif scan.
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
        // overlap=2; annotated at L=94 must beat canonical motif at L=92.
        final ChrIntron annotated = new ChrIntron(CHR1, 1095, 1499);
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        bases[1092] = 'G'; bases[1093] = 'T';
        bases[1495] = 'A'; bases[1496] = 'G';

        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = motifResolver(
                annotated(annotated), java.util.Collections.singletonMap(CHR1, bases)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M405N57M", res.MergedCigar);
        assertEquals(annotated, res.IntroducedIntrons.get(0));
    }

    @Test
    public void testMotifScanPrefersCanonicalOverSemiCanonical()
    {
        // overlap=2; canonical (GT-AG) at L=92 must beat semi-canonical (GC-AG) at L=94.
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        bases[1094] = 'G'; bases[1095] = 'C';   // semi at (1095, 1499): GC-AG
        bases[1497] = 'A'; bases[1498] = 'G';
        bases[1092] = 'G'; bases[1093] = 'T';   // canonical at (1093, 1497): GT-AG
        bases[1495] = 'A'; bases[1496] = 'G';

        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = motifResolver(
                java.util.Collections.emptySet(), java.util.Collections.singletonMap(CHR1, bases)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("92M405N59M", res.MergedCigar);
        assertEquals(new ChrIntron(CHR1, 1093, 1497), res.IntroducedIntrons.get(0));
    }

    @Test
    public void testMotifScanAcceptsReverseStrandCanonical()
    {
        // Reverse-strand transcript: genomic forward CT...AC (RC of GT-AG). Merge via motif scan.
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        bases[1094] = 'C'; bases[1095] = 'T';
        bases[1497] = 'A'; bases[1498] = 'C';

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
        // All-N ref, no annotation, AnnotatedOnly=false — falls through to trust-primary fallback.
        final byte[] bases = new byte[2000];
        java.util.Arrays.fill(bases, (byte) 'N');
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = motifResolver(
                java.util.Collections.emptySet(), java.util.Collections.singletonMap(CHR1, bases)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("94M405N57M", res.MergedCigar);
    }

    @Test
    public void testMotifScanLeftExtendCanonical()
    {
        // Left-extend with canonical GT-AG motif at intron (1058, 1499).
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
        // Trailing 5S "CCCCC" matches chr1:201..205 exactly across annotated intron (131, 200).
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        for(int i = 200; i < 205; ++i) chr1Ref[i] = 'C';

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

        assertTrue(res.Merged);
        assertEquals("30M70N5M", res.MergedCigar);
        assertEquals(101, res.MergedStart);
    }

    @Test
    public void testRefVerifyRightExtendSnapsBackOverExtendedBoundary()
    {
        // bwa over-extended 1 base into the intron (31M4S); boundary snap trims back to find
        // the annotated intron at 131 and ref-verifies the 5-base tail.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        chr1Ref[130] = 'C';
        for(int i = 200; i < 205; ++i) chr1Ref[i] = 'C';

        final byte[] readBases = new byte[35];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';
        for(int i = 30; i < 35; ++i) readBases[i] = 'C';

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 35, 101, "31M4S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("30M70N5M", res.MergedCigar);
        assertEquals(101, res.MergedStart);
    }

    @Test
    public void testRefVerifyBothEndsClippedRescuesJunctionTailKeepsOtherClip()
    {
        // Both ends clipped: trailing 5S is the junction tail; leading 5S must survive untouched.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        for(int i = 200; i < 205; ++i) chr1Ref[i] = 'C';

        final byte[] readBases = new byte[40];
        for(int i = 0; i < 5; ++i) readBases[i] = 'T';
        for(int i = 5; i < 35; ++i) readBases[i] = 'A';
        for(int i = 35; i < 40; ++i) readBases[i] = 'C';

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 40, 101, "5S30M5S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("5S30M70N5M", res.MergedCigar);
        assertEquals(101, res.MergedStart);
    }

    @Test
    public void testRefVerifyRightExtendRejectsOnMismatch()
    {
        // Ref all 'G'; read has 'C' — mismatch rejects the merge.
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
        // 20-base softclip: tolerance floor(20/10)=2; 1 mismatch — accept.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        for(int i = 200; i < 218; ++i) chr1Ref[i] = 'C';
        chr1Ref[218] = 'G';
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
        // Leading 5S "TTTTT" matches chr1:126..130 exactly across annotated intron (131, 200).
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
        // Empty annotation — no candidate exon for the trailing softclip.
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
        // Two annotated introns share donor; both downstream exons match equally — ambiguous.
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
        // No RefSequenceSource — ref-verify skips silently.
        final RescueCandidate cand = candidate(CHR1, true, 35, 101, "30M5S", 60);

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.NO_MATCHING_SUPP, res.RejectReason);
    }

    @Test
    public void testRefVerifyPartialMatchTrailingKeepsOuterClip()
    {
        // 15 proximal bases match the exon; outer 4 are adapter residual — only proximal bases convert
        // to M, outer stay clipped: 30M19S → 30M70N15M4S.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        for(int i = 200; i < 215; ++i) chr1Ref[i] = 'C';    // downstream exon: 15 bases

        final byte[] readBases = new byte[49];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';
        for(int i = 30; i < 45; ++i) readBases[i] = 'C';    // 15-base junction overhang
        for(int i = 45; i < 49; ++i) readBases[i] = 'T';    // 4-base adapter residual

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 49, 101, "30M19S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("30M70N15M4S", res.MergedCigar);
        assertEquals(101, res.MergedStart);
    }

    @Test
    public void testRefVerifyPartialMatchLeadingKeepsOuterClip()
    {
        // Mirror on leading side: outer 4 are adapter, proximal 15 are real overhang → 4S15M70N30M.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        for(int i = 115; i < 130; ++i) chr1Ref[i] = 'C';    // upstream exon: 15 bases

        final byte[] readBases = new byte[49];
        for(int i = 0; i < 4; ++i) readBases[i] = 'T';      // 4-base adapter residual
        for(int i = 4; i < 19; ++i) readBases[i] = 'C';     // 15-base junction overhang
        for(int i = 19; i < 49; ++i) readBases[i] = 'A';

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 49, 201, "19S30M", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("4S15M70N30M", res.MergedCigar);
        assertEquals(116, res.MergedStart);
    }

    @Test
    public void testRefVerifyShortPartialRunRejected()
    {
        // 8-base proximal match < MinPartialMatchRun (11) — guards against chance match manufacturing a junction.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        for(int i = 200; i < 208; ++i) chr1Ref[i] = 'C';    // only 8 matching bases

        final byte[] readBases = new byte[49];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';
        for(int i = 30; i < 38; ++i) readBases[i] = 'C';    // 8-base proximal match
        for(int i = 38; i < 49; ++i) readBases[i] = 'T';    // 11-base divergent residual

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 49, 101, "30M19S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertFalse(res.Merged);
        assertEquals(RescueRejectReason.REF_VERIFY_SHORT_PARTIAL_RUN, res.RejectReason);
    }

    @Test
    public void testRefVerifySnapsBackSevenBaseOverExtension()
    {
        // bwa over-extended 7 bases into the intron (37M3S); boundary snap (MaxBoundaryShift≥8)
        // trims back to find annotated intron 131 and re-verifies the 10-base tail.
        final byte[] chr1Ref = new byte[300];
        for(int i = 0; i < chr1Ref.length; ++i) chr1Ref[i] = 'A';
        for(int i = 130; i < 137; ++i) chr1Ref[i] = 'C';    // bwa's 7 over-extended bases
        for(int i = 200; i < 210; ++i) chr1Ref[i] = 'C';    // real downstream exon

        final byte[] readBases = new byte[40];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';    // 7 over-extended + 3 clipped

        final ChrIntron intron = new ChrIntron(CHR1, 131, 200);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 40, 101, "37M3S", 60, Collections.emptyList(),
                readBases, Collections.emptyList());

        final RescueResult res = refVerifyResolver(
                new java.util.HashSet<>(Arrays.asList(intron)),
                java.util.Collections.singletonMap(CHR1, chr1Ref)).resolve(cand);

        assertTrue(res.Merged);
        assertEquals("30M70N10M", res.MergedCigar);
        assertEquals(101, res.MergedStart);
    }

    @Test
    public void testMateHintUsedAsFallbackWhenAnnotationMisses()
    {
        // Mate hint overrides trust-primary fallback when no annotation is within the snap window.
        final RescueCandidate withoutHint = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));
        final RescueResult resNoHint = enabledResolver(annotated()).resolve(withoutHint);
        assertTrue(resNoHint.Merged);
        assertEquals(1095, resNoHint.IntroducedIntrons.get(0).IntronStart);

        final ChrIntron hint = new ChrIntron(CHR1, 1093, 1497);
        final RescueCandidate withHint = new RescueCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));
        final RescueResult resHinted = enabledResolver(annotated()).resolve(withHint);
        assertTrue(resHinted.Merged);
        assertEquals(1093, resHinted.IntroducedIntrons.get(0).IntronStart);
        assertEquals("92M405N59M", resHinted.MergedCigar);
    }

    @Test
    public void testMateHintIgnoredWhenAnnotatedMatchExists()
    {
        // Annotated junction beats mate hint when both are within the snap window.
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
        // Hint outside the overlap window — ignored, falls back to trust-primary.
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
