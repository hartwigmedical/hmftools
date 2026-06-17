package com.hartwig.hmftools.tars.liftback.rescue;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.repeatedBase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.SpliceCommon;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

// Tests for JunctionRescueResolver. Reads are 151bp (Illumina default) unless noted.
public class JunctionRescueResolverTest
{
    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";
    private static final int READ_LEN = 151;

    private static Set<ChrBaseRegion> annotated(final ChrBaseRegion... introns)
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

    private JunctionRescueResolver enabledResolver(final Set<ChrBaseRegion> annotated)
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

        final ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 31448462, 31448540);
        final RescueSupplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M79N57M", res.mergedCigar());
        assertEquals(primStart, res.mergedStart());
        assertEquals(1, res.droppedSupplementaryIndices().size());
        assertEquals(Integer.valueOf(0), res.droppedSupplementaryIndices().get(0));
        assertEquals(1, res.introducedIntrons().size());
        assertEquals(annotatedIntron, res.introducedIntrons().get(0));
        assertEquals(1, res.chainDepth());
    }

    @Test
    public void testFoldTrailingFabricatedMicroJunction()
    {
        // 100M83N3M48S: 3bp anchor split off by a spurious 83N against the trailing clip -> 100M51S.
        assertEquals("100M51S", CigarUtils.cigarElementsToStr(JunctionRescueResolver.foldFabricatedTerminalMicroJunctions(
                CigarUtils.cigarElementsFromStr("100M83N3M48S"), SpliceCommon.MIN_JUNCTION_ANCHOR)));
    }

    @Test
    public void testFoldLeadingFabricatedMicroJunction()
    {
        assertEquals("51S100M", CigarUtils.cigarElementsToStr(JunctionRescueResolver.foldFabricatedTerminalMicroJunctions(
                CigarUtils.cigarElementsFromStr("48S3M83N100M"), SpliceCommon.MIN_JUNCTION_ANCHOR)));
    }

    @Test
    public void testNoFoldWhenAnchorAboveThreshold()
    {
        // 8bp anchor meets MIN_JUNCTION_ANCHOR (8) - a trusted junction, left intact.
        assertEquals("100M83N8M48S", CigarUtils.cigarElementsToStr(JunctionRescueResolver.foldFabricatedTerminalMicroJunctions(
                CigarUtils.cigarElementsFromStr("100M83N8M48S"), SpliceCommon.MIN_JUNCTION_ANCHOR)));
    }

    @Test
    public void testNoFoldWithoutTerminalSoftClip()
    {
        // no terminal softclip - nothing to fold into.
        assertEquals("100M83N3M", CigarUtils.cigarElementsToStr(JunctionRescueResolver.foldFabricatedTerminalMicroJunctions(
                CigarUtils.cigarElementsFromStr("100M83N3M"), SpliceCommon.MIN_JUNCTION_ANCHOR)));
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

        final ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1051370, 1051525);
        final RescueSupplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("100M156N51M", res.mergedCigar());
        assertEquals(primStart, res.mergedStart());
        assertEquals(annotatedIntron, res.introducedIntrons().get(0));
    }

    @Test
    public void testLeftExtendCleanMerge()
    {
        // Mirror of right-extend: primary starts after junction, supp covers the upstream exon.
        final int suppStart = 31448368;
        final String suppCigar = "57M94S";
        final int primStart = 31448541;
        final String primCigar = "57S94M";

        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 31448425, 31448540);
        final RescueSupplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("57M116N94M", res.mergedCigar());
        assertEquals(suppStart, res.mergedStart());
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

        final Set<ChrBaseRegion> set = annotated(
                new ChrBaseRegion(CHR1, 1050, 1999),
                new ChrBaseRegion(CHR1, 2060, 2999));

        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255,
                supp(0, CHR1, true, suppMidStart, suppMidCigar, 60),
                supp(1, CHR1, true, suppLastStart, suppLastCigar, 60));

        final RescueResult res = enabledResolver(set).resolve(cand);

        assertTrue(res.merged());
        assertEquals("50M950N60M940N41M", res.mergedCigar());
        assertEquals(primStart, res.mergedStart());
        assertEquals(2, res.chainDepth());
        assertEquals(2, res.droppedSupplementaryIndices().size());
        assertTrue(res.droppedSupplementaryIndices().contains(0));
        assertTrue(res.droppedSupplementaryIndices().contains(1));
    }

    @Test
    public void testMergeWhenPrimaryAlreadyHasInternalN()
    {
        // Primary already has an internal junction; supp picks up the next exon.
        final int primStart = 1000;
        final String primCigar = "50M200N40M61S";
        final int suppStart = 1500;
        final String suppCigar = "90S61M";

        final Set<ChrBaseRegion> set = annotated(new ChrBaseRegion(CHR1, 1290, 1499));

        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255,
                supp(0, CHR1, true, suppStart, suppCigar, 60));

        final RescueResult res = enabledResolver(set).resolve(cand);

        assertTrue(res.merged());
        assertEquals("50M200N40M210N61M", res.mergedCigar());
    }

    @Test
    public void testRejectWhenDisabled()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1175, "90S61M", 60));

        final RescueResult res = disabledResolver().resolve(cand);

        assertFalse(res.merged());
    }

    @Test
    public void testRejectNoTerminalSoftclip()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "151M", 60,
                supp(0, CHR1, true, 2000, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.NO_TERMINAL_SOFTCLIP, res.rejectReason());
    }

    @Test
    public void testRejectDifferentChromosome()
    {
        final ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR2, true, 1500, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.DIFFERENT_CHROMOSOME, res.rejectReason());
    }

    @Test
    public void testRejectOppositeStrand()
    {
        final ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, false, 1500, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.OPPOSITE_STRAND, res.rejectReason());
    }

    @Test
    public void testRejectIntronTooShort()
    {
        // Intron length 6 - below default MinIntronLength=21.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1100, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.INTRON_TOO_SHORT, res.rejectReason());
    }

    @Test
    public void testRejectIntronTooLong()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 2_001_095, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.INTRON_TOO_LONG, res.rejectReason());
    }

    @Test
    public void testRejectShortPrimaryAnchor()
    {
        // Primary anchor 2bp - below MinAnchorOverhang=3.
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1002, 1099);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "2M149S", 60,
                supp(0, CHR1, true, 1100, "2S149M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.SHORT_ANCHOR, res.rejectReason());
    }

    @Test
    public void testRejectShortSuppAnchor()
    {
        // Supp anchor 1bp - below MinAnchorOverhang=3.
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1050, 1099);
        final RescueCandidate cand = candidate(CHR1, true, 50, 1000, "49M1S", 60,
                supp(0, CHR1, true, 1100, "49S1M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.SHORT_ANCHOR, res.rejectReason());
    }

    @Test
    public void testRejectNovelJunctionWhenAnnotatedOnly()
    {
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.NOVEL_JUNCTION, res.rejectReason());
    }

    @Test
    public void testAcceptNovelJunctionWhenAnnotatedOnlyFalse()
    {
        final RescueConfig perm = new RescueConfig(true, 21, 1_000_000, 3, 4, false, 0, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), perm).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M106N57M", res.mergedCigar());
    }

    @Test
    public void testRejectHardClipPrimary()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "10H94M57S", 60,
                supp(0, CHR1, true, 1200, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.COMPLEX_CIGAR_SHAPE, res.rejectReason());
    }

    @Test
    public void testRejectHardClipSupp()
    {
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1199);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1200, "5H90S61M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.COMPLEX_CIGAR_SHAPE, res.rejectReason());
    }

    @Test
    public void testRejectIndelAdjacentToPrimarySoftclip()
    {
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1091, 1199);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "90M4I57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.COMPLEX_CIGAR_SHAPE, res.rejectReason());
    }

    @Test
    public void testRejectReadCoverageOverlap()
    {
        // overlap 14 bases - exceeds tolerance.
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "80S71M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.READ_COVERAGE_OVERLAP, res.rejectReason());
    }

    @Test
    public void testRejectReadCoverageGap()
    {
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "110S41M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.READ_COVERAGE_GAP, res.rejectReason());
    }

    @Test
    public void testRejectRefOverlap()
    {
        // Supp starts upstream of primary's matched end - ref overlap.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1080, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.READ_COVERAGE_OVERLAP, res.rejectReason());
    }

    @Test
    public void testNoSupplementaryAvailable()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60);

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.NO_MATCHING_SUPP, res.rejectReason());
    }

    @Test
    public void testRejectShapeMismatchSuppWrongClipSide()
    {
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "61M90S", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.NO_MATCHING_SUPP, res.rejectReason());
    }

    @Test
    public void testPrimaryBothSidesClippedRightExtendAccepted()
    {
        // Middle-anchored primary: supp past primaryRefEnd disambiguates direction → right-extend fires.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                supp(0, CHR1, true, 1500, "95S56M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.merged());
        assertEquals("5S90M409N56M", res.mergedCigar());
        assertEquals(1001, res.mergedStart());
    }

    @Test
    public void testPrimaryBothSidesClippedLeftExtendAccepted()
    {
        // Mirror on the left: supp ends before primaryStart → left-extend fires.
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                supp(0, CHR1, true, 500, "5M146S", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.merged());
        assertEquals("5M496N90M56S", res.mergedCigar());
        assertEquals(500, res.mergedStart());
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

        assertTrue(res.merged());
        assertEquals("5M496N90M409N56M", res.mergedCigar());
        assertEquals(500, res.mergedStart());
        assertEquals(2, res.chainDepth());
    }

    @Test
    public void testMultipleSuppsInReachSkipsMerge()
    {
        // Two supps both pass every gate - refuse to guess the splice destination.
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.MULTIPLE_SUPPS_IN_REACH, res.rejectReason());
    }

    @Test
    public void testPrimaryMapqZeroStillMerges()
    {
        // MAPQ-0 primary is the common tx-contig duplicate artifact - not a disqualifier.
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 31448462, 31448540);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 31448368, "94M57S", 0,
                supp(0, CHR1, true, 31448541, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.merged());
    }

    @Test
    public void testTwoValidSuppsSameBoundaryRefusesToGuess()
    {
        // Two supps at different intron lengths both pass - single-supp-within-reach policy refuses to guess.
        final ChrBaseRegion intron1 = new ChrBaseRegion(CHR1, 1095, 1499);
        final ChrBaseRegion intron2 = new ChrBaseRegion(CHR1, 1095, 1799);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1800, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron1, intron2)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.MULTIPLE_SUPPS_IN_REACH, res.rejectReason());
    }

    @Test
    public void testMergeWhenSuppMapqZero()
    {
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 0));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.merged());
    }

    @Test
    public void testChainDepthCap()
    {
        // cap=2 stops the chain after 2 merges even when more supps are available.
        final RescueConfig cappedConfig = new RescueConfig(true, 21, 1_000_000, 3, 2, true, 0, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);

        final int p = 1000;
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, p, "30M121S", 60,
                supp(0, CHR1, true, 2000, "30S30M91S", 60),
                supp(1, CHR1, true, 3000, "60S30M61S", 60),
                supp(2, CHR1, true, 4000, "90S30M31S", 60),
                supp(3, CHR1, true, 5000, "120S31M", 60));

        final Set<ChrBaseRegion> set = annotated(
                new ChrBaseRegion(CHR1, 1030, 1999),
                new ChrBaseRegion(CHR1, 2030, 2999),
                new ChrBaseRegion(CHR1, 3030, 3999),
                new ChrBaseRegion(CHR1, 4030, 4999));

        final RescueResult res = new JunctionRescueResolver(set, cappedConfig).resolve(cand);

        assertTrue(res.merged());
        assertEquals(2, res.chainDepth());
        assertEquals(2, res.droppedSupplementaryIndices().size());
    }

    @Test
    public void testShortReadLength()
    {
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1031, 1130);
        final RescueCandidate cand = candidate(CHR1, true, 50, 1001, "30M20S", 60,
                supp(0, CHR1, true, 1131, "30S20M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("30M100N20M", res.mergedCigar());
    }

    @Test
    public void testStatisticsCounters()
    {
        // AnnotatedOnly=true so the novel-junction reject is observable.
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final JunctionRescueResolver resolver = new JunctionRescueResolver(
                annotated(new ChrBaseRegion(CHR1, 1095, 1499)), strict);

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

        assertFalse(res.merged());
        assertEquals(1, resolver.statistics().candidatesEvaluated());
        assertEquals(0, resolver.statistics().mergedTotal());
        assertEquals(1, resolver.statistics().rejectCount(RescueRejectReason.OPPOSITE_STRAND));
    }

    @Test
    public void testExp7Case2Chr1_31448368()
    {
        // exp7 case 2 (chr1:31448368): clean complementary cigars merge into the expected junction CIGAR.
        final ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 31448462, 31448539);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 31448368, "94M57S", 60,
                supp(0, CHR1, true, 31448540, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M78N57M", res.mergedCigar());
        assertEquals(31448368, res.mergedStart());
    }

    @Test
    public void testExp7Case3Chr5_34937631()
    {
        // exp7 chr5:34937631: overlap=2, snap picks L=61 (trust-supp) on the annotated junction,
        // producing the expected 61M1166N90M.
        final String chr5 = "chr5";
        final ChrBaseRegion annotatedIntron = new ChrBaseRegion(chr5, 34937692, 34938857);
        final RescueCandidate cand = candidate(chr5, true, READ_LEN, 34938856, "59S92M", 60,
                supp(0, chr5, true, 34937631, "61M90S", 60));

        final RescueResult res = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("61M1166N90M", res.mergedCigar());
        assertEquals(34937631, res.mergedStart());
    }

    @Test
    public void testOverlapWithinToleranceSnapsToAnnotatedJunction()
    {
        // overlap=4; trust-primary L=94 lands on annotated (1095, 1503).
        final ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1095, 1503);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "90S61M", 60));

        final RescueResult res = enabledResolver(annotated(annotated)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M409N57M", res.mergedCigar());
        assertEquals(1001, res.mergedStart());
        assertEquals(annotated, res.introducedIntrons().get(0));
    }

    @Test
    public void testOverlapWithinToleranceSnapsToAnnotatedAtTrustSuppEnd()
    {
        // overlap=2; trust-supp L=92 lands on annotated (1092, 1497).
        final ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1092, 1497);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = enabledResolver(annotated(annotated)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("92M406N59M", res.mergedCigar());
    }

    @Test
    public void testOverlapWithinToleranceMaxMinAnchorWins()
    {
        // overlap=2, both L's annotated. max-min-anchor: trust-supp (min=59) beats trust-primary (min=55).
        final ChrBaseRegion trustPrimary = new ChrBaseRegion(CHR1, 1095, 1499);
        final ChrBaseRegion trustSupp = new ChrBaseRegion(CHR1, 1093, 1497);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = enabledResolver(annotated(trustPrimary, trustSupp)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("92M405N59M", res.mergedCigar());
        assertEquals(trustSupp, res.introducedIntrons().get(0));
    }

    @Test
    public void testOverlapExceedsToleranceRejected()
    {
        final ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "88S63M", 60));

        final RescueResult res = enabledResolver(annotated(annotated)).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.READ_COVERAGE_OVERLAP, res.rejectReason());
    }

    @Test
    public void testNoAnnotatedSnapWithAnnotatedOnlyTrueRejects()
    {
        final RescueConfig strict = new RescueConfig(true, 21, 1_000_000, 3, 4, true, 5, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.NOVEL_JUNCTION, res.rejectReason());
    }

    @Test
    public void testNoAnnotatedSnapWithDefaultsFallsBackToTrustPrimary()
    {
        // AnnotatedOnly=false - merge succeeds via trust-primary fallback when no L is annotated.
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M407N57M", res.mergedCigar());
    }

    @Test
    public void testNoAnnotatedSnapWithAnnotatedOnlyFalseFallsBackToTrustPrimary()
    {
        final RescueConfig perm = new RescueConfig(true, 21, 1_000_000, 3, 4, false, 5, 0, RescueConfig.DEFAULT_MIN_PARTIAL_MATCH_RUN);
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        final RescueResult res = new JunctionRescueResolver(Collections.emptySet(), perm).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M407N57M", res.mergedCigar());
    }

    @Test
    public void testMergedIntronListReflectsChainOrder()
    {
        final int primStart = 1000;
        final Set<ChrBaseRegion> set = annotated(
                new ChrBaseRegion(CHR1, 1050, 1999),
                new ChrBaseRegion(CHR1, 2060, 2999));

        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, primStart, "50M101S", 60,
                supp(0, CHR1, true, 2000, "50S60M41S", 60),
                supp(1, CHR1, true, 3000, "110S41M", 60));

        final RescueResult res = enabledResolver(set).resolve(cand);

        assertTrue(res.merged());
        assertEquals(2, res.introducedIntrons().size());
        assertEquals(new ChrBaseRegion(CHR1, 1050, 1999), res.introducedIntrons().get(0));
        assertEquals(new ChrBaseRegion(CHR1, 2060, 2999), res.introducedIntrons().get(1));
    }

    // Resolver wired to a base-level genome, for the motif-scan and ref-verify passes.
    private static JunctionRescueResolver resolverWithRef(final Set<ChrBaseRegion> annotated, final TestGenome genome)
    {
        return new JunctionRescueResolver(
                new AnnotatedJunctionIndex(annotated), genome.asRefSource(), RescueConfig.enabledDefaults());
    }

    // 'N' genome with a canonical GT-AG motif seeded at the intron flanks.
    private static TestGenome refWithCanonicalIntron(final int chromLen, final int intronStart, final int intronEnd)
    {
        return new TestGenome().with(CHR1, chromLen, 'N')
                .set(CHR1, intronStart, "GT").set(CHR1, intronEnd - 1, "AG");
    }

    // Candidate carrying read bases for ref-verify: no supps, no mate hints.
    private static RescueCandidate refVerifyCandidate(final int start, final String cigar, final int readLen, final byte[] readBases)
    {
        return new RescueCandidate(CHR1, true, readLen, start, cigar, 60,
                Collections.emptyList(), readBases, Collections.emptyList());
    }

    @Test
    public void testMotifScanPicksCanonicalGTagWhenUnannotated()
    {
        // Canonical GT-AG motif placed at intron (1095, 1499) with no annotation - merge via motif scan.
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = resolverWithRef(
                Collections.emptySet(), refWithCanonicalIntron(2000, 1095, 1499)).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M405N57M", res.mergedCigar());
        assertEquals(new ChrBaseRegion(CHR1, 1095, 1499), res.introducedIntrons().get(0));
    }

    @Test
    public void testMotifScanPicksAnnotatedOverMotifWhenBothPresent()
    {
        // overlap=2; annotated at L=94 must beat canonical motif at L=92.
        final ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1095, 1499);
        final TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1093, "GT").set(CHR1, 1496, "AG");

        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = resolverWithRef(annotated(annotated), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M405N57M", res.mergedCigar());
        assertEquals(annotated, res.introducedIntrons().get(0));
    }

    @Test
    public void testMotifScanPrefersCanonicalOverSemiCanonical()
    {
        // overlap=2; canonical (GT-AG) at L=92 must beat semi-canonical (GC-AG) at L=94.
        final TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1095, "GC").set(CHR1, 1498, "AG")    // semi at (1095, 1499): GC-AG
                .set(CHR1, 1093, "GT").set(CHR1, 1496, "AG");   // canonical at (1093, 1497): GT-AG

        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        final RescueResult res = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("92M405N59M", res.mergedCigar());
        assertEquals(new ChrBaseRegion(CHR1, 1093, 1497), res.introducedIntrons().get(0));
    }

    @Test
    public void testMotifScanAcceptsReverseStrandCanonical()
    {
        // Reverse-strand transcript: genomic forward CT...AC (RC of GT-AG). Merge via motif scan.
        final TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1095, "CT").set(CHR1, 1498, "AC");

        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M405N57M", res.mergedCigar());
    }

    @Test
    public void testMotifScanIgnoredWhenNoMotifAndNoAnnotated()
    {
        // All-N ref, no annotation, AnnotatedOnly=false - falls through to trust-primary fallback.
        final TestGenome genome = new TestGenome().with(CHR1, 2000, 'N');
        final RescueCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("94M405N57M", res.mergedCigar());
    }

    @Test
    public void testMotifScanLeftExtendCanonical()
    {
        // Left-extend with canonical GT-AG motif at intron (1058, 1499).
        final TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1058, "GT").set(CHR1, 1498, "AG");

        final RescueCandidate cand = candidate(CHR1, true, 151, 1500, "57S94M", 60,
                supp(0, CHR1, true, 1001, "57M94S", 60));

        final RescueResult res = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("57M442N94M", res.mergedCigar());
    }

    @Test
    public void testRefVerifyRightExtendSuccess()
    {
        // Trailing 5S "CCCCC" matches chr1:201..205 exactly across annotated intron (131, 200).
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 5, 'C');
        final byte[] readBases = bases("A".repeat(30) + "C".repeat(5));

        final RescueCandidate cand = refVerifyCandidate(101, "30M5S", 35, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("30M70N5M", res.mergedCigar());
        assertEquals(101, res.mergedStart());
    }

    @Test
    public void testRefVerifyRightExtendSnapsBackOverExtendedBoundary()
    {
        // bwa over-extended 1 base into the intron (31M4S); boundary snap trims back to find
        // the annotated intron at 131 and ref-verifies the 5-base tail.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A')
                .set(CHR1, 131, "C").set(CHR1, 201, 5, 'C');
        final byte[] readBases = bases("A".repeat(30) + "C".repeat(5));

        final RescueCandidate cand = refVerifyCandidate(101, "31M4S", 35, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("30M70N5M", res.mergedCigar());
        assertEquals(101, res.mergedStart());
    }

    @Test
    public void testRefVerifyBothEndsClippedRescuesJunctionTailKeepsOtherClip()
    {
        // Both ends clipped: trailing 5S is the junction tail; leading 5S must survive untouched.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 5, 'C');
        final byte[] readBases = bases("T".repeat(5) + "A".repeat(30) + "C".repeat(5));

        final RescueCandidate cand = refVerifyCandidate(101, "5S30M5S", 40, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("5S30M70N5M", res.mergedCigar());
        assertEquals(101, res.mergedStart());
    }

    @Test
    public void testRefVerifyRightExtendRejectsOnMismatch()
    {
        // Ref all 'G'; read has 'C' - mismatch rejects the merge.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'G');
        final byte[] readBases = bases("A".repeat(30) + "C".repeat(5));

        final RescueCandidate cand = refVerifyCandidate(101, "30M5S", 35, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH, res.rejectReason());
    }

    @Test
    public void testRefVerifyClipsMismatchNotRecoveredByScore()
    {
        // 20-base trailing clip: 18 proximal matches, then a mismatch with only 1 trailing match. Under the
        // shared bwa-mem score walk (mismatch -4) one trailing match cannot recover the mismatch, so the run
        // stops at 18 and the outer 2 bases stay soft-clipped.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A')
                .set(CHR1, 201, 18, 'C').set(CHR1, 219, "GC");
        final byte[] readBases = bases("A".repeat(20) + "C".repeat(20));

        final RescueCandidate cand = refVerifyCandidate(101, "20M20S", 40, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 121, 200)), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("20M80N18M2S", res.mergedCigar());
    }

    @Test
    public void testRefVerifyLeftExtendSuccess()
    {
        // Leading 5S "TTTTT" matches chr1:126..130 exactly across annotated intron (131, 200).
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 126, 5, 'T');
        final byte[] readBases = bases("T".repeat(5) + "A".repeat(30));

        final RescueCandidate cand = refVerifyCandidate(201, "5S30M", 35, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("5M70N30M", res.mergedCigar());
        assertEquals(126, res.mergedStart());
    }

    @Test
    public void testRefVerifyRejectsWhenNoCandidateExon()
    {
        // Empty annotation - no candidate exon for the trailing softclip.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A');
        final RescueCandidate cand = refVerifyCandidate(101, "30M5S", 35, repeatedBase(35, 'A'));

        final RescueResult res = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.REF_VERIFY_NO_CANDIDATE_EXON, res.rejectReason());
    }

    @Test
    public void testRefVerifyAmbiguousRejected()
    {
        // Two annotated introns share donor; both downstream exons (201..205 and 501..505) match "CCCCC" - ambiguous.
        final TestGenome genome = new TestGenome().with(CHR1, 600, 'A')
                .set(CHR1, 201, 5, 'C').set(CHR1, 501, 5, 'C');
        final byte[] readBases = bases("A".repeat(30) + "C".repeat(5));

        final RescueCandidate cand = refVerifyCandidate(101, "30M5S", 35, readBases);

        final RescueResult res = resolverWithRef(
                annotated(new ChrBaseRegion(CHR1, 131, 200), new ChrBaseRegion(CHR1, 131, 500)), genome).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.REF_VERIFY_AMBIGUOUS, res.rejectReason());
    }

    @Test
    public void testRefVerifyWithoutRefSourceSkipsSilently()
    {
        // No RefSequenceSource - ref-verify skips silently.
        final RescueCandidate cand = candidate(CHR1, true, 35, 101, "30M5S", 60);

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.NO_MATCHING_SUPP, res.rejectReason());
    }

    @Test
    public void testRefVerifyPartialMatchTrailingKeepsOuterClip()
    {
        // 15 proximal bases match the exon; outer 4 are adapter residual - only proximal bases convert
        // to M, outer stay clipped: 30M19S → 30M70N15M4S.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 15, 'C');   // downstream exon: 15 bases
        final byte[] readBases = bases("A".repeat(30) + "C".repeat(15) + "T".repeat(4));   // overhang + adapter residual

        final RescueCandidate cand = refVerifyCandidate(101, "30M19S", 49, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("30M70N15M4S", res.mergedCigar());
        assertEquals(101, res.mergedStart());
    }

    @Test
    public void testRefVerifyPartialMatchLeadingKeepsOuterClip()
    {
        // Mirror on leading side: outer 4 are adapter, proximal 15 are real overhang → 4S15M70N30M.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 116, 15, 'C');   // upstream exon: 15 bases
        final byte[] readBases = bases("T".repeat(4) + "C".repeat(15) + "A".repeat(30));   // adapter + overhang

        final RescueCandidate cand = refVerifyCandidate(201, "19S30M", 49, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("4S15M70N30M", res.mergedCigar());
        assertEquals(116, res.mergedStart());
    }

    @Test
    public void testRefVerifyShortPartialRunRejected()
    {
        // 8-base proximal match < MinPartialMatchRun (11) - guards against chance match manufacturing a junction.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 8, 'C');   // only 8 matching bases
        final byte[] readBases = bases("A".repeat(30) + "C".repeat(8) + "T".repeat(11));   // 8-base match + divergent residual

        final RescueCandidate cand = refVerifyCandidate(101, "30M19S", 49, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertFalse(res.merged());
        assertEquals(RescueRejectReason.REF_VERIFY_SHORT_PARTIAL_RUN, res.rejectReason());
    }

    @Test
    public void testRefVerifySnapsBackSevenBaseOverExtension()
    {
        // bwa over-extended 7 bases into the intron (37M3S); boundary snap (MaxBoundaryShift≥8)
        // trims back to find annotated intron 131 and re-verifies the 10-base tail.
        final TestGenome genome = new TestGenome().with(CHR1, 300, 'A')
                .set(CHR1, 131, 7, 'C')     // bwa's 7 over-extended bases
                .set(CHR1, 201, 10, 'C');   // real downstream exon
        final byte[] readBases = bases("A".repeat(30) + "C".repeat(10));   // 7 over-extended + 3 clipped

        final RescueCandidate cand = refVerifyCandidate(101, "37M3S", 40, readBases);

        final RescueResult res = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(res.merged());
        assertEquals("30M70N10M", res.mergedCigar());
        assertEquals(101, res.mergedStart());
    }

    @Test
    public void testMateHintUsedAsFallbackWhenAnnotationMisses()
    {
        // Mate hint overrides trust-primary fallback when no annotation is within the snap window.
        final RescueCandidate withoutHint = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));
        final RescueResult resNoHint = enabledResolver(annotated()).resolve(withoutHint);
        assertTrue(resNoHint.merged());
        assertEquals(1095, resNoHint.introducedIntrons().get(0).start());

        final ChrBaseRegion hint = new ChrBaseRegion(CHR1, 1093, 1497);
        final RescueCandidate withHint = new RescueCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));
        final RescueResult resHinted = enabledResolver(annotated()).resolve(withHint);
        assertTrue(resHinted.merged());
        assertEquals(1093, resHinted.introducedIntrons().get(0).start());
        assertEquals("92M405N59M", resHinted.mergedCigar());
    }

    @Test
    public void testMateHintIgnoredWhenAnnotatedMatchExists()
    {
        // Annotated junction beats mate hint when both are within the snap window.
        final ChrBaseRegion annotatedJunc = new ChrBaseRegion(CHR1, 1095, 1499);
        final ChrBaseRegion hint = new ChrBaseRegion(CHR1, 1093, 1497);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));

        final RescueResult res = enabledResolver(annotated(annotatedJunc)).resolve(cand);

        assertTrue(res.merged());
        assertEquals(1095, res.introducedIntrons().get(0).start());
    }

    @Test
    public void testMateHintOutsideOverlapWindowIgnored()
    {
        // Hint outside the overlap window - ignored, falls back to trust-primary.
        final ChrBaseRegion hint = new ChrBaseRegion(CHR1, 50000, 50100);
        final RescueCandidate cand = new RescueCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));

        final RescueResult res = enabledResolver(annotated()).resolve(cand);

        assertTrue(res.merged());
        assertEquals(1095, res.introducedIntrons().get(0).start());
    }

    @Test
    public void testRejectReasonIsNullOnSuccess()
    {
        final ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        final RescueCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        final RescueResult res = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(res.merged());
        assertNull(res.rejectReason());
    }
}
