package com.hartwig.hmftools.tars.liftback.supplementary;

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

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

// Tests for SupplementaryResolver. Reads are 151bp (Illumina default) unless noted.
public class SupplementaryResolverTest
{
    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";
    private static final int READ_LEN = 151;

    private static Set<ChrBaseRegion> annotated(final ChrBaseRegion... introns)
    {
        return new HashSet<>(Arrays.asList(introns));
    }

    private static SupplementaryCandidate candidate(
            final String chrom, final boolean forward, final int readLen, final int primStart,
            final String primCigar, final int primMapq, final SupplementaryRecord... supps)
    {
        return new SupplementaryCandidate(chrom, forward, readLen, primStart, primCigar, primMapq,
                supps.length == 0 ? Collections.emptyList() : Arrays.asList(supps));
    }

    private static SupplementaryRecord supp(
            final int index, final String chrom, final boolean forward, final int start,
            final String cigar, final int mapq)
    {
        return new SupplementaryRecord(index, chrom, forward, start, cigar, mapq);
    }

    private SupplementaryResolver disabledResolver()
    {
        return new SupplementaryResolver(Collections.emptySet(), SupplementaryConfig.defaults());
    }

    private SupplementaryResolver enabledResolver(final Set<ChrBaseRegion> annotated)
    {
        return new SupplementaryResolver(annotated, SupplementaryConfig.enabledDefaults());
    }

    @Test
    public void testRightExtendCleanMerge()
    {
        // Clean complementary cigars across an annotated intron; tests core right-extend.
        int primStart = 31448368;
        String primCigar = "94M57S";
        int suppStart = 31448541;
        String suppCigar = "94S57M";

        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 31448462, 31448540);
        SupplementaryRecord supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        SupplementaryResult result = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M79N57M", result.mergedCigar());
        assertEquals(primStart, result.mergedStart());
        assertEquals(1, result.droppedSupplementaryIndices().size());
        assertEquals(Integer.valueOf(0), result.droppedSupplementaryIndices().get(0));
        assertEquals(1, result.introducedIntrons().size());
        assertEquals(annotatedIntron, result.introducedIntrons().get(0));
        assertEquals(1, result.chainDepth());
    }

    @Test
    public void testResolveAcrossPreCollapsedTerminalClip()
    {
        // exp8 read 25535: the overhang gate now collapses the tx-contig over-run (100M83N3M48S) to 100M51S
        // upstream, so supplementary resolve receives the clean clip and merges it across the true 156N junction.
        // (The fold that used to do this inside the resolver was removed; the gate owns terminal micro-junction handling.)
        int primStart = 1051270;
        String primCigar = "100M51S";
        int suppStart = 1051525;
        String suppCigar = "99S52M";

        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1051370, 1051525);
        SupplementaryRecord supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        SupplementaryResult result = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("100M156N51M", result.mergedCigar());
        assertEquals(primStart, result.mergedStart());
        assertEquals(annotatedIntron, result.introducedIntrons().get(0));
    }

    @Test
    public void testLeftExtendCleanMerge()
    {
        // Mirror of right-extend: primary starts after junction, supp covers the upstream exon.
        int suppStart = 31448368;
        String suppCigar = "57M94S";
        int primStart = 31448541;
        String primCigar = "57S94M";

        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 31448425, 31448540);
        SupplementaryRecord supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255, supp);

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("57M116N94M", result.mergedCigar());
        assertEquals(suppStart, result.mergedStart());
    }

    @Test
    public void testChainMergeAcrossThreeExons()
    {
        // 3-exon read: primary on first exon, two supps on middle and last exons.
        int primStart = 1000;
        String primCigar = "50M101S";
        int suppMidStart = 2000;
        String suppMidCigar = "50S60M41S";
        int suppLastStart = 3000;
        String suppLastCigar = "110S41M";

        Set<ChrBaseRegion> set = annotated(
                new ChrBaseRegion(CHR1, 1050, 1999),
                new ChrBaseRegion(CHR1, 2060, 2999));

        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255,
                supp(0, CHR1, true, suppMidStart, suppMidCigar, 60),
                supp(1, CHR1, true, suppLastStart, suppLastCigar, 60));

        SupplementaryResult result = enabledResolver(set).resolve(cand);

        assertTrue(result.merged());
        assertEquals("50M950N60M940N41M", result.mergedCigar());
        assertEquals(primStart, result.mergedStart());
        assertEquals(2, result.chainDepth());
        assertEquals(2, result.droppedSupplementaryIndices().size());
        assertTrue(result.droppedSupplementaryIndices().contains(0));
        assertTrue(result.droppedSupplementaryIndices().contains(1));

        // introduced introns are recorded in chain order (first exon boundary first).
        assertEquals(2, result.introducedIntrons().size());
        assertEquals(new ChrBaseRegion(CHR1, 1050, 1999), result.introducedIntrons().get(0));
        assertEquals(new ChrBaseRegion(CHR1, 2060, 2999), result.introducedIntrons().get(1));
    }

    @Test
    public void testMergeWhenPrimaryAlreadyHasInternalN()
    {
        // Primary already has an internal junction; supp picks up the next exon.
        int primStart = 1000;
        String primCigar = "50M200N40M61S";
        int suppStart = 1500;
        String suppCigar = "90S61M";

        Set<ChrBaseRegion> set = annotated(new ChrBaseRegion(CHR1, 1290, 1499));

        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, 255,
                supp(0, CHR1, true, suppStart, suppCigar, 60));

        SupplementaryResult result = enabledResolver(set).resolve(cand);

        assertTrue(result.merged());
        assertEquals("50M200N40M210N61M", result.mergedCigar());
    }

    @Test
    public void testRejectWhenDisabled()
    {
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1175, "90S61M", 60));

        SupplementaryResult result = disabledResolver().resolve(cand);

        assertFalse(result.merged());
    }

    @Test
    public void testRejectNoTerminalSoftclip()
    {
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "151M", 60,
                supp(0, CHR1, true, 2000, "90S61M", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.NO_TERMINAL_SOFTCLIP, result.rejectReason());
    }

    @Test
    public void testRejectDifferentChromosome()
    {
        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1095, 1499);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR2, true, 1500, "90S61M", 60));

        SupplementaryResult result = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.DIFFERENT_CHROMOSOME, result.rejectReason());
    }

    @Test
    public void testRejectOppositeStrand()
    {
        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1095, 1499);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, false, 1500, "90S61M", 60));

        SupplementaryResult result = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.OPPOSITE_STRAND, result.rejectReason());
    }

    @Test
    public void testRejectIntronTooShort()
    {
        // Intron length 6 - below default MinIntronLength=21.
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1100, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.INTRON_TOO_SHORT, result.rejectReason());
    }

    @Test
    public void testRejectIntronTooLong()
    {
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 2_001_095, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.INTRON_TOO_LONG, result.rejectReason());
    }

    @Test
    public void testShortPrimaryAnchorNowMerges()
    {
        // Primary anchor 2bp: the MinAnchorOverhang guard was removed, so this merges on the annotated junction.
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1002, 1099);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "2M149S", 60,
                supp(0, CHR1, true, 1100, "2S149M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("2M98N149M", result.mergedCigar());
    }

    @Test
    public void testShortSuppAnchorNowMerges()
    {
        // Supp anchor 1bp: the MinAnchorOverhang guard was removed, so this merges (trust-primary junction).
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1050, 1099);
        SupplementaryCandidate cand = candidate(CHR1, true, 50, 1000, "49M1S", 60,
                supp(0, CHR1, true, 1100, "49S1M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("49M51N1M", result.mergedCigar());
    }

    @Test
    public void testRejectNovelJunctionWhenAnnotatedOnly()
    {
        SupplementaryConfig strict = new SupplementaryConfig(true, 21, 1_000_000, 4, true, 5, 0);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        SupplementaryResult result = new SupplementaryResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.NOVEL_JUNCTION, result.rejectReason());
    }

    @Test
    public void testAcceptNovelJunctionWhenAnnotatedOnlyFalse()
    {
        SupplementaryConfig perm = new SupplementaryConfig(true, 21, 1_000_000, 4, false, 0, 0);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        SupplementaryResult result = new SupplementaryResolver(Collections.emptySet(), perm).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M106N57M", result.mergedCigar());
    }

    @Test
    public void testRejectHardClipPrimary()
    {
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "10H94M57S", 60,
                supp(0, CHR1, true, 1200, "90S61M", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE, result.rejectReason());
    }

    @Test
    public void testRejectHardClipSupp()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1199);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1200, "5H90S61M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE, result.rejectReason());
    }

    @Test
    public void testRejectIndelAdjacentToPrimarySoftclip()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1091, 1199);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "90M4I57S", 60,
                supp(0, CHR1, true, 1200, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.COMPLEX_CIGAR_SHAPE, result.rejectReason());
    }

    @Test
    public void testRejectReadCoverageOverlap()
    {
        // overlap 14 bases - exceeds tolerance.
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "80S71M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.READ_COVERAGE_OVERLAP, result.rejectReason());
    }

    @Test
    public void testRejectReadCoverageGap()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "110S41M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.READ_COVERAGE_GAP, result.rejectReason());
    }

    @Test
    public void testRejectRefOverlap()
    {
        // Supp starts upstream of primary's matched end - ref overlap.
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1080, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.READ_COVERAGE_OVERLAP, result.rejectReason());
    }

    @Test
    public void testNoSupplementaryAvailable()
    {
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60);

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.NO_MATCHING_SUPP, result.rejectReason());
    }

    @Test
    public void testRejectShapeMismatchSuppWrongClipSide()
    {
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1200, "61M90S", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.NO_MATCHING_SUPP, result.rejectReason());
    }

    @Test
    public void testPrimaryBothSidesClippedRightExtendAccepted()
    {
        // Middle-anchored primary: supp past primaryRefEnd disambiguates direction -> right-extend fires.
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                supp(0, CHR1, true, 1500, "95S56M", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals("5S90M409N56M", result.mergedCigar());
        assertEquals(1001, result.mergedStart());
    }

    @Test
    public void testPrimaryBothSidesClippedLeftExtendAccepted()
    {
        // Mirror on the left: supp ends before primaryStart -> left-extend fires.
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                supp(0, CHR1, true, 500, "5M146S", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals("5M496N90M56S", result.mergedCigar());
        assertEquals(500, result.mergedStart());
    }

    @Test
    public void testPrimaryBothSidesClippedChainMergesBothSupps()
    {
        // Full middle-anchored 3-exon scenario: chain merges right first then left.
        SupplementaryRecord right = supp(0, CHR1, true, 1500, "95S56M", 60);
        SupplementaryRecord left = supp(1, CHR1, true, 500, "5M146S", 60);
        SupplementaryCandidate cand = new SupplementaryCandidate(CHR1, true, READ_LEN, 1001, "5S90M56S", 60,
                Arrays.asList(right, left));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals("5M496N90M409N56M", result.mergedCigar());
        assertEquals(500, result.mergedStart());
        assertEquals(2, result.chainDepth());
    }

    @Test
    public void testMultipleSuppsInReachSkipsMerge()
    {
        // Two supps both pass every gate - refuse to guess the splice destination.
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1500, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.MULTIPLE_SUPPS_IN_REACH, result.rejectReason());
    }

    @Test
    public void testPrimaryMapqZeroStillMerges()
    {
        // MAPQ-0 primary is the common tx-contig duplicate artifact - not a disqualifier.
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 31448462, 31448540);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 31448368, "94M57S", 0,
                supp(0, CHR1, true, 31448541, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
    }

    @Test
    public void testTwoValidSuppsSameBoundaryRefusesToGuess()
    {
        // Two supps at different intron lengths both pass - single-supp-within-reach policy refuses to guess.
        ChrBaseRegion intron1 = new ChrBaseRegion(CHR1, 1095, 1499);
        ChrBaseRegion intron2 = new ChrBaseRegion(CHR1, 1095, 1799);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1800, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron1, intron2)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.MULTIPLE_SUPPS_IN_REACH, result.rejectReason());
    }

    @Test
    public void testMergeWhenSuppMapqZero()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 0));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
    }

    @Test
    public void testChainDepthCap()
    {
        // cap=2 stops the chain after 2 merges even when more supps are available.
        SupplementaryConfig cappedConfig =
                new SupplementaryConfig(true, 21, 1_000_000, 2, true, 0, 0);

        int primaryStart = 1000;
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, primaryStart, "30M121S", 60,
                supp(0, CHR1, true, 2000, "30S30M91S", 60),
                supp(1, CHR1, true, 3000, "60S30M61S", 60),
                supp(2, CHR1, true, 4000, "90S30M31S", 60),
                supp(3, CHR1, true, 5000, "120S31M", 60));

        Set<ChrBaseRegion> set = annotated(
                new ChrBaseRegion(CHR1, 1030, 1999),
                new ChrBaseRegion(CHR1, 2030, 2999),
                new ChrBaseRegion(CHR1, 3030, 3999),
                new ChrBaseRegion(CHR1, 4030, 4999));

        SupplementaryResult result = new SupplementaryResolver(set, cappedConfig).resolve(cand);

        assertTrue(result.merged());
        assertEquals(2, result.chainDepth());
        assertEquals(2, result.droppedSupplementaryIndices().size());
    }

    @Test
    public void testShortReadLength()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1031, 1130);
        SupplementaryCandidate cand = candidate(CHR1, true, 50, 1001, "30M20S", 60,
                supp(0, CHR1, true, 1131, "30S20M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("30M100N20M", result.mergedCigar());
    }

    @Test
    public void testStatisticsCounters()
    {
        // AnnotatedOnly=true so the novel-junction reject is observable.
        SupplementaryConfig strict = new SupplementaryConfig(true, 21, 1_000_000, 4, true, 5, 0);
        SupplementaryResolver resolver = new SupplementaryResolver(
                annotated(new ChrBaseRegion(CHR1, 1095, 1499)), strict);

        resolver.resolve(candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60)));
        resolver.resolve(candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1800, "94S57M", 60)));
        resolver.resolve(candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR2, true, 1500, "94S57M", 60)));

        SupplementaryStatistics stats = resolver.statistics();
        assertEquals(3, stats.candidatesEvaluated());
        assertEquals(1, stats.mergedTotal());
        assertEquals(1, stats.mergedAtChainDepth(1));
        assertEquals(1, stats.rejectCount(SupplementaryRejectReason.NOVEL_JUNCTION));
        assertEquals(1, stats.rejectCount(SupplementaryRejectReason.DIFFERENT_CHROMOSOME));
    }

    @Test
    public void testNoMergeStillCountsCandidateAndReject()
    {
        SupplementaryResolver resolver = new SupplementaryResolver(
                Collections.emptySet(), SupplementaryConfig.enabledDefaults());

        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S", 60,
                supp(0, CHR1, false, 1500, "94S57M", 60));

        SupplementaryResult result = resolver.resolve(cand);

        assertFalse(result.merged());
        assertEquals(1, resolver.statistics().candidatesEvaluated());
        assertEquals(0, resolver.statistics().mergedTotal());
        assertEquals(1, resolver.statistics().rejectCount(SupplementaryRejectReason.OPPOSITE_STRAND));
    }

    @Test
    public void testExp7Case2Chr1_31448368()
    {
        // exp7 case 2 (chr1:31448368): clean complementary cigars merge into the expected junction CIGAR.
        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 31448462, 31448539);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 31448368, "94M57S", 60,
                supp(0, CHR1, true, 31448540, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M78N57M", result.mergedCigar());
        assertEquals(31448368, result.mergedStart());
    }

    @Test
    public void testExp7Case3Chr5_34937631()
    {
        // exp7 chr5:34937631: overlap=2, snap picks L=61 (trust-supp) on the annotated junction,
        // producing the expected 61M1166N90M.
        String chr5 = "chr5";
        ChrBaseRegion annotatedIntron = new ChrBaseRegion(chr5, 34937692, 34938857);
        SupplementaryCandidate cand = candidate(chr5, true, READ_LEN, 34938856, "59S92M", 60,
                supp(0, chr5, true, 34937631, "61M90S", 60));

        SupplementaryResult result = enabledResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("61M1166N90M", result.mergedCigar());
        assertEquals(34937631, result.mergedStart());
    }

    @Test
    public void testOverlapWithinToleranceSnapsToAnnotatedJunction()
    {
        // overlap=4; trust-primary L=94 lands on annotated (1095, 1503).
        ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1095, 1503);
        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "90S61M", 60));

        SupplementaryResult result = enabledResolver(annotated(annotated)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M409N57M", result.mergedCigar());
        assertEquals(1001, result.mergedStart());
        assertEquals(annotated, result.introducedIntrons().get(0));
    }

    @Test
    public void testOverlapWithinToleranceSnapsToAnnotatedAtTrustSuppEnd()
    {
        // overlap=2; trust-supp L=92 lands on annotated (1092, 1497).
        ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1092, 1497);
        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1000, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        SupplementaryResult result = enabledResolver(annotated(annotated)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("92M406N59M", result.mergedCigar());
    }

    @Test
    public void testOverlapWithinToleranceSeededPickAmongAnnotated()
    {
        // overlap=2, both L's land on an annotated junction (a tie at the ANNOTATED tier). The pick is seeded by
        // the read (deterministic), not by max-min-anchor; this candidate resolves to trustPrimary.
        ChrBaseRegion trustPrimary = new ChrBaseRegion(CHR1, 1095, 1499);
        ChrBaseRegion trustSupp = new ChrBaseRegion(CHR1, 1093, 1497);
        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        SupplementaryResult result = enabledResolver(annotated(trustPrimary, trustSupp)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M405N57M", result.mergedCigar());
        assertEquals(trustPrimary, result.introducedIntrons().get(0));
    }

    @Test
    public void testOverlapExceedsToleranceRejected()
    {
        ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1095, 1499);
        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "88S63M", 60));

        SupplementaryResult result = enabledResolver(annotated(annotated)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.READ_COVERAGE_OVERLAP, result.rejectReason());
    }

    @Test
    public void testNoAnnotatedSnapWithAnnotatedOnlyTrueRejects()
    {
        SupplementaryConfig strict = new SupplementaryConfig(true, 21, 1_000_000, 4, true, 5, 0);
        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        SupplementaryResult result = new SupplementaryResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.NOVEL_JUNCTION, result.rejectReason());
    }

    @Test
    public void testNoAnnotatedSnapWithDefaultsFallsBackToMidpoint()
    {
        // AnnotatedOnly=false - with no annotated/motif snap the junction is placed at the midpoint of the
        // ambiguous overlap range (rounded down), not at bwa's split point.
        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "92S59M", 60));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals("93M407N58M", result.mergedCigar());
    }

    // Resolver wired to a base-level genome, for the motif-scan and ref-verify passes.
    private static SupplementaryResolver resolverWithRef(final Set<ChrBaseRegion> annotated, final TestGenome genome)
    {
        return new SupplementaryResolver(
                new AnnotatedJunctionIndex(annotated), genome.asRefSource(), SupplementaryConfig.enabledDefaults());
    }

    // 'N' genome with a canonical GT-AG motif seeded at the intron flanks.
    private static TestGenome refWithCanonicalIntron(final int chromLen, final int intronStart, final int intronEnd)
    {
        return new TestGenome().with(CHR1, chromLen, 'N')
                .set(CHR1, intronStart, "GT").set(CHR1, intronEnd - 1, "AG");
    }

    // Candidate carrying read bases for ref-verify: no supps, no mate hints.
    private static SupplementaryCandidate refVerifyCandidate(final int start, final String cigar, final int readLen, final byte[] readBases)
    {
        return new SupplementaryCandidate(CHR1, true, readLen, start, cigar, 60,
                Collections.emptyList(), readBases, Collections.emptyList());
    }

    @Test
    public void testMotifScanPicksCanonicalGTagWhenUnannotated()
    {
        // Canonical GT-AG motif placed at intron (1095, 1499) with no annotation - merge via motif scan.
        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        SupplementaryResult result = resolverWithRef(
                Collections.emptySet(), refWithCanonicalIntron(2000, 1095, 1499)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M405N57M", result.mergedCigar());
        assertEquals(new ChrBaseRegion(CHR1, 1095, 1499), result.introducedIntrons().get(0));
    }

    @Test
    public void testMotifScanPicksAnnotatedOverMotifWhenBothPresent()
    {
        // overlap=2; annotated at L=94 must beat canonical motif at L=92.
        ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1095, 1499);
        TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1093, "GT").set(CHR1, 1496, "AG");

        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        SupplementaryResult result = resolverWithRef(annotated(annotated), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M405N57M", result.mergedCigar());
        assertEquals(annotated, result.introducedIntrons().get(0));
    }

    @Test
    public void testMotifScanPrefersCanonicalOverSemiCanonical()
    {
        // overlap=2; canonical (GT-AG) at L=92 must beat semi-canonical (GC-AG) at L=94.
        TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1095, "GC").set(CHR1, 1498, "AG")    // semi at (1095, 1499): GC-AG
                .set(CHR1, 1093, "GT").set(CHR1, 1496, "AG");   // canonical at (1093, 1497): GT-AG

        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));

        SupplementaryResult result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("92M405N59M", result.mergedCigar());
        assertEquals(new ChrBaseRegion(CHR1, 1093, 1497), result.introducedIntrons().get(0));
    }

    @Test
    public void testMotifScanAcceptsReverseStrandCanonical()
    {
        // Reverse-strand transcript: genomic forward CT...AC (RC of GT-AG). Merge via motif scan.
        TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1095, "CT").set(CHR1, 1498, "AC");

        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        SupplementaryResult result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M405N57M", result.mergedCigar());
    }

    @Test
    public void testMotifScanIgnoredWhenNoMotifAndNoAnnotated()
    {
        // All-N ref, no annotation, AnnotatedOnly=false - falls through to trust-primary fallback.
        TestGenome genome = new TestGenome().with(CHR1, 2000, 'N');
        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        SupplementaryResult result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M405N57M", result.mergedCigar());
    }

    @Test
    public void testMotifScanLeftExtendCanonical()
    {
        // Left-extend with canonical GT-AG motif at intron (1058, 1499).
        TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1058, "GT").set(CHR1, 1498, "AG");

        SupplementaryCandidate cand = candidate(CHR1, true, 151, 1500, "57S94M", 60,
                supp(0, CHR1, true, 1001, "57M94S", 60));

        SupplementaryResult result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("57M442N94M", result.mergedCigar());
    }

    @Test
    public void testRefVerifyRightExtendSuccess()
    {
        // Trailing 5S "CCCCC" matches chr1:201..205 exactly across annotated intron (131, 200).
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 5, 'C');
        byte[] readBases = bases("A".repeat(30) + "C".repeat(5));

        SupplementaryCandidate cand = refVerifyCandidate(101, "30M5S", 35, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("30M70N5M", result.mergedCigar());
        assertEquals(101, result.mergedStart());
    }

    @Test
    public void testRefVerifyRightExtendSnapsBackOverExtendedBoundary()
    {
        // bwa over-extended 1 base into the intron (31M4S); boundary snap trims back to find
        // the annotated intron at 131 and ref-verifies the 5-base tail.
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A')
                .set(CHR1, 131, "C").set(CHR1, 201, 5, 'C');
        byte[] readBases = bases("A".repeat(30) + "C".repeat(5));

        SupplementaryCandidate cand = refVerifyCandidate(101, "31M4S", 35, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("30M70N5M", result.mergedCigar());
        assertEquals(101, result.mergedStart());
    }

    @Test
    public void testRefVerifyBothEndsClippedResolvesJunctionTailKeepsOtherClip()
    {
        // Both ends clipped: trailing 5S is the junction tail; leading 5S must survive untouched.
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 5, 'C');
        byte[] readBases = bases("T".repeat(5) + "A".repeat(30) + "C".repeat(5));

        SupplementaryCandidate cand = refVerifyCandidate(101, "5S30M5S", 40, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("5S30M70N5M", result.mergedCigar());
        assertEquals(101, result.mergedStart());
    }

    @Test
    public void testRefVerifyRightExtendRejectsOnMismatch()
    {
        // Ref all 'G'; read has 'C' - mismatch rejects the merge.
        TestGenome genome = new TestGenome().with(CHR1, 300, 'G');
        byte[] readBases = bases("A".repeat(30) + "C".repeat(5));

        SupplementaryCandidate cand = refVerifyCandidate(101, "30M5S", 35, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.REF_VERIFY_MISMATCH_TOO_HIGH, result.rejectReason());
    }

    @Test
    public void testRefVerifyClipsMismatchNotRecoveredByScore()
    {
        // 20-base trailing clip: 18 proximal matches, then a mismatch with only 1 trailing match. Under the
        // shared bwa-mem score walk (mismatch -4) one trailing match cannot recover the mismatch, so the run
        // stops at 18 and the outer 2 bases stay soft-clipped.
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A')
                .set(CHR1, 201, 18, 'C').set(CHR1, 219, "GC");
        byte[] readBases = bases("A".repeat(20) + "C".repeat(20));

        SupplementaryCandidate cand = refVerifyCandidate(101, "20M20S", 40, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 121, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("20M80N18M2S", result.mergedCigar());
    }

    @Test
    public void testRefVerifyLeftExtendSuccess()
    {
        // Leading 5S "TTTTT" matches chr1:126..130 exactly across annotated intron (131, 200).
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 126, 5, 'T');
        byte[] readBases = bases("T".repeat(5) + "A".repeat(30));

        SupplementaryCandidate cand = refVerifyCandidate(201, "5S30M", 35, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("5M70N30M", result.mergedCigar());
        assertEquals(126, result.mergedStart());
    }

    @Test
    public void testRefVerifyRejectsWhenNoCandidateExon()
    {
        // Empty annotation - no candidate exon for the trailing softclip.
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A');
        SupplementaryCandidate cand = refVerifyCandidate(101, "30M5S", 35, repeatedBase(35, 'A'));

        SupplementaryResult result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.REF_VERIFY_NO_CANDIDATE_EXON, result.rejectReason());
    }

    @Test
    public void testRefVerifyAmbiguousRejected()
    {
        // Two annotated introns share donor; both downstream exons (201..205 and 501..505) match "CCCCC" - ambiguous.
        TestGenome genome = new TestGenome().with(CHR1, 600, 'A')
                .set(CHR1, 201, 5, 'C').set(CHR1, 501, 5, 'C');
        byte[] readBases = bases("A".repeat(30) + "C".repeat(5));

        SupplementaryCandidate cand = refVerifyCandidate(101, "30M5S", 35, readBases);

        SupplementaryResult result = resolverWithRef(
                annotated(new ChrBaseRegion(CHR1, 131, 200), new ChrBaseRegion(CHR1, 131, 500)), genome).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.REF_VERIFY_AMBIGUOUS, result.rejectReason());
    }

    @Test
    public void testRefVerifyWithoutRefSourceSkipsSilently()
    {
        // No RefSequenceSource - ref-verify skips silently.
        SupplementaryCandidate cand = candidate(CHR1, true, 35, 101, "30M5S", 60);

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(SupplementaryRejectReason.NO_MATCHING_SUPP, result.rejectReason());
    }

    @Test
    public void testRefVerifyPartialMatchTrailingKeepsOuterClip()
    {
        // 15 proximal bases match the exon; outer 4 are adapter residual - only proximal bases convert
        // to M, outer stay clipped: 30M19S -> 30M70N15M4S.
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 15, 'C');   // downstream exon: 15 bases
        byte[] readBases = bases("A".repeat(30) + "C".repeat(15) + "T".repeat(4));   // overhang + adapter residual

        SupplementaryCandidate cand = refVerifyCandidate(101, "30M19S", 49, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("30M70N15M4S", result.mergedCigar());
        assertEquals(101, result.mergedStart());
    }

    @Test
    public void testRefVerifyPartialMatchLeadingKeepsOuterClip()
    {
        // Mirror on leading side: outer 4 are adapter, proximal 15 are real overhang -> 4S15M70N30M.
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 116, 15, 'C');   // upstream exon: 15 bases
        byte[] readBases = bases("T".repeat(4) + "C".repeat(15) + "A".repeat(30));   // adapter + overhang

        SupplementaryCandidate cand = refVerifyCandidate(201, "19S30M", 49, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("4S15M70N30M", result.mergedCigar());
        assertEquals(116, result.mergedStart());
    }

    @Test
    public void testRefVerifyShortPartialRunNowMerges()
    {
        // 8-base proximal match: the MinPartialMatchRun guard was removed, so the partial match now merges with
        // the divergent residual left soft-clipped (30M70N8M11S).
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 8, 'C');   // only 8 matching bases
        byte[] readBases = bases("A".repeat(30) + "C".repeat(8) + "T".repeat(11));   // 8-base match + divergent residual

        SupplementaryCandidate cand = refVerifyCandidate(101, "30M19S", 49, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("30M70N8M11S", result.mergedCigar());
        assertEquals(101, result.mergedStart());
    }

    @Test
    public void testRefVerifySnapsBackSevenBaseOverExtension()
    {
        // bwa over-extended 7 bases into the intron (37M3S); boundary snap (MaxBoundaryShift>=8)
        // trims back to find annotated intron 131 and re-verifies the 10-base tail.
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A')
                .set(CHR1, 131, 7, 'C')     // bwa's 7 over-extended bases
                .set(CHR1, 201, 10, 'C');   // real downstream exon
        byte[] readBases = bases("A".repeat(30) + "C".repeat(10));   // 7 over-extended + 3 clipped

        SupplementaryCandidate cand = refVerifyCandidate(101, "37M3S", 40, readBases);

        SupplementaryResult result = resolverWithRef(annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("30M70N10M", result.mergedCigar());
        assertEquals(101, result.mergedStart());
    }

    @Test
    public void testMateHintUsedAsFallbackWhenAnnotationMisses()
    {
        // Mate hint overrides the midpoint fallback when no annotation is within the snap window.
        SupplementaryCandidate withoutHint = candidate(CHR1, true, 151, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1498, "92S59M", 60));
        SupplementaryResult resNoHint = enabledResolver(annotated()).resolve(withoutHint);
        assertTrue(resNoHint.merged());
        assertEquals(1094, resNoHint.introducedIntrons().get(0).start());

        ChrBaseRegion hint = new ChrBaseRegion(CHR1, 1093, 1497);
        SupplementaryCandidate withHint = new SupplementaryCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));
        SupplementaryResult resHinted = enabledResolver(annotated()).resolve(withHint);
        assertTrue(resHinted.merged());
        assertEquals(1093, resHinted.introducedIntrons().get(0).start());
        assertEquals("92M405N59M", resHinted.mergedCigar());
    }

    @Test
    public void testMateHintIgnoredWhenAnnotatedMatchExists()
    {
        // Annotated junction beats mate hint when both are within the snap window.
        ChrBaseRegion annotatedJunc = new ChrBaseRegion(CHR1, 1095, 1499);
        ChrBaseRegion hint = new ChrBaseRegion(CHR1, 1093, 1497);
        SupplementaryCandidate cand = new SupplementaryCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));

        SupplementaryResult result = enabledResolver(annotated(annotatedJunc)).resolve(cand);

        assertTrue(result.merged());
        assertEquals(1095, result.introducedIntrons().get(0).start());
    }

    @Test
    public void testMateHintOutsideOverlapWindowIgnored()
    {
        // Hint outside the overlap window - ignored, falls back to the midpoint of the range.
        ChrBaseRegion hint = new ChrBaseRegion(CHR1, 50000, 50100);
        SupplementaryCandidate cand = new SupplementaryCandidate(
                CHR1, true, 151, 1001, "94M57S", 60,
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));

        SupplementaryResult result = enabledResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals(1094, result.introducedIntrons().get(0).start());
    }

    @Test
    public void testRejectReasonIsNullOnSuccess()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        SupplementaryCandidate cand = candidate(CHR1, true, READ_LEN, 1001, "94M57S", 60,
                supp(0, CHR1, true, 1500, "94S57M", 60));

        SupplementaryResult result = enabledResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertNull(result.rejectReason());
    }
}
