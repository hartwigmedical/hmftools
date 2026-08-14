package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.repeatedBase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.liftback.EnsemblAnnotationIndex;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver.Candidate;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver.RejectReason;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver.Result;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver.Supplementary;
import com.hartwig.hmftools.tars.liftback.features.SupplementaryResolver.Tier;

import org.junit.Test;

// Reads are 151bp (Illumina default) unless noted.
public class SupplementaryResolverTest
{
    private static final String CHR1 = "chr1";
    private static final String CHR2 = "chr2";
    private static final int READ_LEN = 151;

    private static Set<ChrBaseRegion> annotated(final ChrBaseRegion... introns)
    {
        return new HashSet<>(Arrays.asList(introns));
    }

    private static Candidate candidate(
            final String chrom, final boolean forward, final int readLen, final int primStart,
            final String primCigar, final Supplementary... supps)
    {
        return new Candidate(
                chrom, forward, readLen, primStart, primCigar,
                supps.length == 0 ? Collections.emptyList() : Arrays.asList(supps));
    }

    private static Supplementary supp(
            final int index, final String chrom, final boolean forward, final int start,
            final String cigar, final int mapQuality)
    {
        return new Supplementary(index, chrom, forward, start, cigar, mapQuality);
    }

    private SupplementaryResolver defaultResolver(final Set<ChrBaseRegion> annotated)
    {
        return new SupplementaryResolver(annotated, SupplementaryConfig.defaults());
    }

    @Test
    public void testRightExtendCleanMerge()
    {
        int primStart = 31448368;
        String primCigar = "94M57S";
        int suppStart = 31448541;
        String suppCigar = "94S57M";

        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 31448462, 31448540);
        Supplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        Candidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, supp);

        Result result = defaultResolver(annotated(annotatedIntron)).resolve(cand);

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
    public void testSpliceStrandComesFromAnnotationAndMotif()
    {
        EnsemblAnnotationIndex annotation = EnsemblAnnotationIndex.fromContigEntries(List.of(
                ContigEntry.annotationOnly(
                        "g", "gn", "tn", CHR1, -1,
                        List.of(new BaseRegion(100, 199), new BaseRegion(300, 399)))));
        SupplementaryResolver annotatedResolver = new SupplementaryResolver(
                annotation, null, SupplementaryConfig.defaults());

        assertEquals(-1, annotatedResolver.spliceStrand(CHR1, 150, "50M100N50M"));
        assertEquals(1, SupplementaryResolver.motifStrand(bases("GT"), bases("AG")));
        assertEquals(-1, SupplementaryResolver.motifStrand(bases("CT"), bases("AC")));
    }

    @Test
    public void testDoesNotInferAJunctionWithoutASupplementary()
    {
        TestGenome genome = new TestGenome().with(CHR1, 300, 'A').set(CHR1, 201, 5, 'C');
        Candidate candidate = new Candidate(
                CHR1, true, 35, 101, "30M5S", Collections.emptyList(),
                bases("A".repeat(30) + "C".repeat(5)), Collections.emptyList());

        Result result = resolverWithRef(
                annotated(new ChrBaseRegion(CHR1, 131, 200)), genome).resolve(candidate);

        assertFalse(result.merged());
        assertEquals(RejectReason.NO_MATCHING_SUPP, result.rejectReason());
    }

    @Test
    public void testResolveAcrossPreCollapsedTerminalClip()
    {
        // The overhang gate collapses the tx-contig over-run (100M83N3M48S) to 100M51S before the resolver runs, so the
        // resolver sees a clean terminal clip and merges it across the true 156N junction.
        int primStart = 1051270;
        String primCigar = "100M51S";
        int suppStart = 1051525;
        String suppCigar = "99S52M";

        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1051370, 1051525);
        Supplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        Candidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, supp);

        Result result = defaultResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("100M156N51M", result.mergedCigar());
        assertEquals(primStart, result.mergedStart());
        assertEquals(annotatedIntron, result.introducedIntrons().get(0));
    }

    @Test
    public void testLeftExtendCleanMerge()
    {
        // Mirror of right-extend: the primary starts after the junction and the supp covers the upstream exon.
        int suppStart = 31448368;
        String suppCigar = "57M94S";
        int primStart = 31448541;
        String primCigar = "57S94M";

        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 31448425, 31448540);
        Supplementary supp = supp(0, CHR1, true, suppStart, suppCigar, 60);
        Candidate cand = candidate(CHR1, true, READ_LEN, primStart, primCigar, supp);

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("57M116N94M", result.mergedCigar());
        assertEquals(suppStart, result.mergedStart());
    }

    @Test
    public void testChainMergeAcrossThreeExons()
    {
        // 3-exon read: primary on the first exon, two supps on the middle and last exons.
        int primStart = 1000;
        String primCigar = "50M101S";
        int suppMidStart = 2000;
        String suppMidCigar = "50S60M41S";
        int suppLastStart = 3000;
        String suppLastCigar = "110S41M";

        Set<ChrBaseRegion> set = annotated(
                new ChrBaseRegion(CHR1, 1050, 1999),
                new ChrBaseRegion(CHR1, 2060, 2999));

        Candidate cand = candidate(
                CHR1, true, READ_LEN, primStart, primCigar,
                supp(0, CHR1, true, suppMidStart, suppMidCigar, 60),
                supp(1, CHR1, true, suppLastStart, suppLastCigar, 60));

        Result result = defaultResolver(set).resolve(cand);

        assertTrue(result.merged());
        assertEquals("50M950N60M940N41M", result.mergedCigar());
        assertEquals(primStart, result.mergedStart());
        assertEquals(2, result.chainDepth());
        assertEquals(2, result.droppedSupplementaryIndices().size());
        assertTrue(result.droppedSupplementaryIndices().contains(0));
        assertTrue(result.droppedSupplementaryIndices().contains(1));

        // Introduced introns are recorded in chain order, first exon boundary first.
        assertEquals(2, result.introducedIntrons().size());
        assertEquals(new ChrBaseRegion(CHR1, 1050, 1999), result.introducedIntrons().get(0));
        assertEquals(new ChrBaseRegion(CHR1, 2060, 2999), result.introducedIntrons().get(1));
    }

    @Test
    public void testMergeWhenPrimaryAlreadyHasInternalN()
    {
        int primStart = 1000;
        String primCigar = "50M200N40M61S";
        int suppStart = 1500;
        String suppCigar = "90S61M";

        Set<ChrBaseRegion> set = annotated(new ChrBaseRegion(CHR1, 1290, 1499));

        Candidate cand = candidate(
                CHR1, true, READ_LEN, primStart, primCigar,
                supp(0, CHR1, true, suppStart, suppCigar, 60));

        Result result = defaultResolver(set).resolve(cand);

        assertTrue(result.merged());
        assertEquals("50M200N40M210N61M", result.mergedCigar());
    }

    @Test
    public void testRejectNoTerminalSoftclip()
    {
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "151M",
                supp(0, CHR1, true, 2000, "90S61M", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.NO_TERMINAL_SOFTCLIP, result.rejectReason());
    }

    @Test
    public void testRejectDifferentChromosome()
    {
        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1095, 1499);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR2, true, 1500, "90S61M", 60));

        Result result = defaultResolver(annotated(annotatedIntron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.DIFFERENT_CHROMOSOME, result.rejectReason());
    }

    @Test
    public void testRejectOppositeStrand()
    {
        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 1095, 1499);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, false, 1500, "90S61M", 60));

        Result result = defaultResolver(annotated(annotatedIntron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.OPPOSITE_STRAND, result.rejectReason());
    }

    @Test
    public void testRejectIntronTooShort()
    {
        // Intron length 6 - below default MinIntronLength=21.
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "94M57S",
                supp(0, CHR1, true, 1100, "94S57M", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.INTRON_TOO_SHORT, result.rejectReason());
    }

    @Test
    public void testRejectIntronTooLong()
    {
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "94M57S",
                supp(0, CHR1, true, 2_001_095, "94S57M", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.INTRON_TOO_LONG, result.rejectReason());
    }

    @Test
    public void testShortPrimaryAnchorNowMerges()
    {
        // Primary anchor 2bp: the MinAnchorOverhang guard was removed, so this merges on the annotated junction.
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1002, 1099);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "2M149S",
                supp(0, CHR1, true, 1100, "2S149M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("2M98N149M", result.mergedCigar());
    }

    @Test
    public void testShortSuppAnchorNowMerges()
    {
        // Supp anchor 1bp: the MinAnchorOverhang guard was removed, so this merges (trust-primary junction).
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1050, 1099);
        Candidate cand = candidate(
                CHR1, true, 50, 1000, "49M1S",
                supp(0, CHR1, true, 1100, "49S1M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("49M51N1M", result.mergedCigar());
    }

    @Test
    public void testRejectNovelJunctionWhenAnnotatedOnly()
    {
        SupplementaryConfig strict = new SupplementaryConfig(21, 1_000_000, 4, true, 5);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "94M57S",
                supp(0, CHR1, true, 1200, "94S57M", 60));

        Result result = new SupplementaryResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.NOVEL_JUNCTION, result.rejectReason());
    }

    @Test
    public void testAcceptNovelJunctionWhenAnnotatedOnlyFalse()
    {
        SupplementaryConfig perm = new SupplementaryConfig(21, 1_000_000, 4, false, 0);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "94M57S",
                supp(0, CHR1, true, 1200, "94S57M", 60));

        Result result = new SupplementaryResolver(Collections.emptySet(), perm).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M106N57M", result.mergedCigar());
    }

    @Test
    public void testRejectHardClipPrimary()
    {
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "10H94M57S",
                supp(0, CHR1, true, 1200, "90S61M", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.COMPLEX_CIGAR_SHAPE, result.rejectReason());
    }

    @Test
    public void testRejectHardClipSupp()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1199);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1200, "5H90S61M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.COMPLEX_CIGAR_SHAPE, result.rejectReason());
    }

    @Test
    public void testRejectIndelAdjacentToPrimarySoftclip()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1091, 1199);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "90M4I57S",
                supp(0, CHR1, true, 1200, "94S57M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.COMPLEX_CIGAR_SHAPE, result.rejectReason());
    }

    @Test
    public void testRejectReadCoverageOverlap()
    {
        // overlap 14 bases - exceeds tolerance.
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "80S71M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.READ_COVERAGE_OVERLAP, result.rejectReason());
    }

    @Test
    public void testRejectReadCoverageGap()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "110S41M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.READ_COVERAGE_GAP, result.rejectReason());
    }

    @Test
    public void testRejectRefOverlap()
    {
        // Supp starts upstream of primary's matched end - ref overlap.
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1080, "94S57M", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.READ_COVERAGE_OVERLAP, result.rejectReason());
    }

    @Test
    public void testNoSupplementaryAvailable()
    {
        Candidate cand = candidate(CHR1, true, READ_LEN, 1000, "94M57S");

        Result result = defaultResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.NO_MATCHING_SUPP, result.rejectReason());
    }

    @Test
    public void testRejectShapeMismatchSuppWrongClipSide()
    {
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "94M57S",
                supp(0, CHR1, true, 1200, "61M90S", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.NO_MATCHING_SUPP, result.rejectReason());
    }

    @Test
    public void testPrimaryBothSidesClippedRightExtendAccepted()
    {
        // Middle-anchored primary: supp past primaryRefEnd disambiguates direction -> right-extend fires.
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "5S90M56S",
                supp(0, CHR1, true, 1500, "95S56M", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals("5S90M409N56M", result.mergedCigar());
        assertEquals(1001, result.mergedStart());
    }

    @Test
    public void testPrimaryBothSidesClippedLeftExtendAccepted()
    {
        // Mirror on the left: supp ends before primaryStart -> left-extend fires.
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "5S90M56S",
                supp(0, CHR1, true, 500, "5M146S", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals("5M496N90M56S", result.mergedCigar());
        assertEquals(500, result.mergedStart());
    }

    @Test
    public void testPrimaryBothSidesClippedChainMergesBothSupps()
    {
        // Middle-anchored 3-exon read: the chain merges right first, then left.
        Supplementary right = supp(0, CHR1, true, 1500, "95S56M", 60);
        Supplementary left = supp(1, CHR1, true, 500, "5M146S", 60);
        Candidate cand = new Candidate(
                CHR1, true, READ_LEN, 1001, "5S90M56S",
                Arrays.asList(right, left));

        Result result = defaultResolver(annotated()).resolve(cand);

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
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1500, "94S57M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.MULTIPLE_SUPPS_IN_REACH, result.rejectReason());
    }

    @Test
    public void testPrimaryMapQualityZeroStillMerges()
    {
        // MAPQ-0 primary is the common tx-contig duplicate artifact - not a disqualifier.
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 31448462, 31448540);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 31448368, "94M57S",
                supp(0, CHR1, true, 31448541, "94S57M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
    }

    @Test
    public void testTwoValidSuppsSameBoundaryRefusesToGuess()
    {
        // Two supps at different intron lengths both pass - single-supp-within-reach policy refuses to guess.
        ChrBaseRegion intron1 = new ChrBaseRegion(CHR1, 1095, 1499);
        ChrBaseRegion intron2 = new ChrBaseRegion(CHR1, 1095, 1799);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "94S57M", 60),
                supp(1, CHR1, true, 1800, "94S57M", 60));

        Result result = defaultResolver(annotated(intron1, intron2)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.MULTIPLE_SUPPS_IN_REACH, result.rejectReason());
    }

    @Test
    public void testMergeWhenSuppMapQualityZero()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "94S57M", 0));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
    }

    @Test
    public void testChainDepthCap()
    {
        // cap=2 stops the chain after 2 merges even when more supps are available.
        SupplementaryConfig cappedConfig =
                new SupplementaryConfig(21, 1_000_000, 2, true, 0);

        int primaryStart = 1000;
        Candidate cand = candidate(
                CHR1, true, READ_LEN, primaryStart, "30M121S",
                supp(0, CHR1, true, 2000, "30S30M91S", 60),
                supp(1, CHR1, true, 3000, "60S30M61S", 60),
                supp(2, CHR1, true, 4000, "90S30M31S", 60),
                supp(3, CHR1, true, 5000, "120S31M", 60));

        Set<ChrBaseRegion> set = annotated(
                new ChrBaseRegion(CHR1, 1030, 1999),
                new ChrBaseRegion(CHR1, 2030, 2999),
                new ChrBaseRegion(CHR1, 3030, 3999),
                new ChrBaseRegion(CHR1, 4030, 4999));

        Result result = new SupplementaryResolver(set, cappedConfig).resolve(cand);

        assertTrue(result.merged());
        assertEquals(2, result.chainDepth());
        assertEquals(2, result.droppedSupplementaryIndices().size());
    }

    @Test
    public void testShortReadLength()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1031, 1130);
        Candidate cand = candidate(
                CHR1, true, 50, 1001, "30M20S",
                supp(0, CHR1, true, 1131, "30S20M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("30M100N20M", result.mergedCigar());
    }

    @Test
    public void testEachOutcomeReportsItsRejectReason()
    {
        // AnnotatedOnly=true so the novel-junction reject is observable.
        SupplementaryConfig strict = new SupplementaryConfig(21, 1_000_000, 4, true, 5);
        SupplementaryResolver resolver = new SupplementaryResolver(
                annotated(new ChrBaseRegion(CHR1, 1095, 1499)), strict);

        Result onAnnotated = resolver.resolve(candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "94S57M", 60)));
        Result novelJunction = resolver.resolve(candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1800, "94S57M", 60)));
        Result crossChromosome = resolver.resolve(candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR2, true, 1500, "94S57M", 60)));

        assertTrue(onAnnotated.merged());
        assertFalse(novelJunction.merged());
        assertEquals(RejectReason.NOVEL_JUNCTION, novelJunction.rejectReason());
        assertFalse(crossChromosome.merged());
        assertEquals(RejectReason.DIFFERENT_CHROMOSOME, crossChromosome.rejectReason());
    }

    @Test
    public void testOppositeStrandSuppRejected()
    {
        SupplementaryResolver resolver = new SupplementaryResolver(
                Collections.emptySet(), SupplementaryConfig.defaults());

        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1000, "94M57S",
                supp(0, CHR1, false, 1500, "94S57M", 60));

        Result result = resolver.resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.OPPOSITE_STRAND, result.rejectReason());
    }

    @Test
    public void testExp7Case2Chr1_31448368()
    {
        // exp7 case 2 (chr1:31448368): clean complementary cigars merge into the expected junction CIGAR.
        ChrBaseRegion annotatedIntron = new ChrBaseRegion(CHR1, 31448462, 31448539);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 31448368, "94M57S",
                supp(0, CHR1, true, 31448540, "94S57M", 60));

        Result result = defaultResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M78N57M", result.mergedCigar());
        assertEquals(31448368, result.mergedStart());
    }

    @Test
    public void testExp7Case3Chr5_34937631()
    {
        // overlap=2: the junction lands at read offset 61 (trust-supp) on the annotated junction, giving 61M1166N90M.
        String chr5 = "chr5";
        ChrBaseRegion annotatedIntron = new ChrBaseRegion(chr5, 34937692, 34938857);
        Candidate cand = candidate(
                chr5, true, READ_LEN, 34938856, "59S92M",
                supp(0, chr5, true, 34937631, "61M90S", 60));

        Result result = defaultResolver(annotated(annotatedIntron)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("61M1166N90M", result.mergedCigar());
        assertEquals(34937631, result.mergedStart());
    }

    @Test
    public void testOverlapWithinTolerancePlacesJunctionOnAnnotated()
    {
        // overlap=4; trust-primary L=94 lands on annotated (1095, 1503).
        ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1095, 1503);
        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "90S61M", 60));

        Result result = defaultResolver(annotated(annotated)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M409N57M", result.mergedCigar());
        assertEquals(1001, result.mergedStart());
        assertEquals(annotated, result.introducedIntrons().get(0));
    }

    @Test
    public void testOverlapWithinTolerancePlacesJunctionOnAnnotatedAtTrustSuppEnd()
    {
        // overlap=2; trust-supp L=92 lands on annotated (1092, 1497).
        ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1092, 1497);
        Candidate cand = candidate(
                CHR1, true, 151, 1000, "94M57S",
                supp(0, CHR1, true, 1498, "92S59M", 60));

        Result result = defaultResolver(annotated(annotated)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("92M406N59M", result.mergedCigar());
    }

    @Test
    public void testOverlapWithinToleranceSeededPickAmongAnnotated()
    {
        // overlap=2 and both L's land on an annotated junction. The ANNOTATED-tier tie is broken by a read-seeded
        // deterministic pick, not by max-min-anchor; this candidate resolves to trustPrimary.
        ChrBaseRegion trustPrimary = new ChrBaseRegion(CHR1, 1095, 1499);
        ChrBaseRegion trustSupp = new ChrBaseRegion(CHR1, 1093, 1497);
        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1498, "92S59M", 60));

        Result result = defaultResolver(annotated(trustPrimary, trustSupp)).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M405N57M", result.mergedCigar());
        assertEquals(trustPrimary, result.introducedIntrons().get(0));
    }

    @Test
    public void testOverlapExceedsToleranceRejected()
    {
        ChrBaseRegion annotated = new ChrBaseRegion(CHR1, 1095, 1499);
        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "88S63M", 60));

        Result result = defaultResolver(annotated(annotated)).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.READ_COVERAGE_OVERLAP, result.rejectReason());
    }

    @Test
    public void testNoAnnotatedPositionWithAnnotatedOnlyTrueRejects()
    {
        SupplementaryConfig strict = new SupplementaryConfig(21, 1_000_000, 4, true, 5);
        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "92S59M", 60));

        Result result = new SupplementaryResolver(Collections.emptySet(), strict).resolve(cand);

        assertFalse(result.merged());
        assertEquals(RejectReason.NOVEL_JUNCTION, result.rejectReason());
    }

    @Test
    public void testNoAnnotatedPositionWithDefaultsFallsBackToMidpoint()
    {
        // With no annotated or motif match the junction is placed at the midpoint of the ambiguous overlap range
        // (rounded down), not at bwa's split point.
        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "92S59M", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals("93M407N58M", result.mergedCigar());
    }

    private static SupplementaryResolver resolverWithRef(final Set<ChrBaseRegion> annotated, final TestGenome genome)
    {
        return new SupplementaryResolver(
                EnsemblAnnotationIndex.fromJunctions(annotated), genome.asRefGenome(), SupplementaryConfig.defaults());
    }

    private static TestGenome refWithCanonicalIntron(final int chromLen, final int intronStart, final int intronEnd)
    {
        return new TestGenome().with(CHR1, chromLen, 'N')
                .set(CHR1, intronStart, "GT").set(CHR1, intronEnd - 1, "AG");
    }

    @Test
    public void testMotifScanPicksCanonicalGTagWhenUnannotated()
    {
        // Canonical GT-AG motif at intron (1095, 1499) with no annotation: merge via motif scan.
        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "94S57M", 60));

        Result result = resolverWithRef(
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

        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1498, "92S59M", 60));

        Result result = resolverWithRef(annotated(annotated), genome).resolve(cand);

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

        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1498, "92S59M", 60));

        Result result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("92M405N59M", result.mergedCigar());
        assertEquals(new ChrBaseRegion(CHR1, 1093, 1497), result.introducedIntrons().get(0));
    }

    @Test
    public void testMotifScanAcceptsReverseStrandCanonical()
    {
        // Reverse-strand transcript: genomic forward CT-AC is the reverse complement of GT-AG.
        TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1095, "CT").set(CHR1, 1498, "AC");

        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "94S57M", 60));

        Result result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M405N57M", result.mergedCigar());
    }

    @Test
    public void testMotifScanIgnoredWhenNoMotifAndNoAnnotated()
    {
        // All-N ref with no annotation falls through to the trust-primary fallback.
        TestGenome genome = new TestGenome().with(CHR1, 2000, 'N');
        Candidate cand = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "94S57M", 60));

        Result result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M405N57M", result.mergedCigar());
    }

    @Test
    public void testMotifScanLeftExtendCanonical()
    {
        // Left-extend with a canonical GT-AG motif at intron (1058, 1499).
        TestGenome genome = new TestGenome().with(CHR1, 2000, 'N')
                .set(CHR1, 1058, "GT").set(CHR1, 1498, "AG");

        Candidate cand = candidate(
                CHR1, true, 151, 1500, "57S94M",
                supp(0, CHR1, true, 1001, "57M94S", 60));

        Result result = resolverWithRef(Collections.emptySet(), genome).resolve(cand);

        assertTrue(result.merged());
        assertEquals("57M442N94M", result.mergedCigar());
    }

    @Test
    public void testSuppWithInternalNBeforePrimaryClampedToItsFirstBlock()
    {
        // ContigTranslator can expand a cross-exon M into M-N-M. This supp already spans the junction and its post-N
        // block lands on the primary, so it carries no soft clip and pairs with nothing. Clamping it back to its first
        // block (57M94S at 1000) makes it the upstream anchor and the merge re-derives the splice.
        Candidate cand = candidate(
                CHR1, true, 151, 2000, "57S94M",
                supp(0, CHR1, true, 1000, "57M943N94M", 60));

        Result result = defaultResolver(annotated(new ChrBaseRegion(CHR1, 1057, 1999))).resolve(cand);

        assertTrue(result.merged());
        assertEquals("57M943N94M", result.mergedCigar());
        assertEquals(1000, result.mergedStart());
        assertEquals(new ChrBaseRegion(CHR1, 1057, 1999), result.introducedIntrons().get(0));
    }

    @Test
    public void testSuppWithInternalNPastPrimaryClampedToItsLastBlock()
    {
        // Mirror: the supp overruns the primary's end rather than starting before it, so it is clamped to its last
        // block (94S57M at 2037) and becomes the downstream anchor.
        Candidate cand = candidate(
                CHR1, true, 151, 1000, "94M57S",
                supp(0, CHR1, true, 1000, "94M943N57M", 60));

        Result result = defaultResolver(annotated(new ChrBaseRegion(CHR1, 1094, 2036))).resolve(cand);

        assertTrue(result.merged());
        assertEquals("94M943N57M", result.mergedCigar());
        assertEquals(1000, result.mergedStart());
        assertEquals(new ChrBaseRegion(CHR1, 1094, 2036), result.introducedIntrons().get(0));
    }

    @Test
    public void testSuppWithInternalNInsidePrimarySpanNotClamped()
    {
        // Nothing to clamp: the supp sits inside the primary's span, so it keeps its M-N-M shape, has no terminal soft
        // clip and finds no merge partner.
        Candidate cand = candidate(
                CHR1, true, 151, 900, "151M",
                supp(0, CHR1, true, 1000, "57M20N94M", 60));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertFalse(result.merged());
    }

    @Test
    public void testMateHintUsedAsFallbackWhenAnnotationMisses()
    {
        // Mate hint overrides the midpoint fallback when no annotation is within the ambiguous window.
        Candidate withoutHint = candidate(
                CHR1, true, 151, 1001, "94M57S",
                supp(0, CHR1, true, 1498, "92S59M", 60));
        Result resNoHint = defaultResolver(annotated()).resolve(withoutHint);
        assertTrue(resNoHint.merged());
        assertEquals(1094, resNoHint.introducedIntrons().get(0).start());

        ChrBaseRegion hint = new ChrBaseRegion(CHR1, 1093, 1497);
        Candidate withHint = new Candidate(
                CHR1, true, 151, 1001, "94M57S",
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));
        Result resHinted = defaultResolver(annotated()).resolve(withHint);
        assertTrue(resHinted.merged());
        assertEquals(1093, resHinted.introducedIntrons().get(0).start());
        assertEquals("92M405N59M", resHinted.mergedCigar());
    }

    @Test
    public void testMateHintIgnoredWhenAnnotatedMatchExists()
    {
        // Annotated junction beats mate hint when both are within the ambiguous window.
        ChrBaseRegion annotatedJunc = new ChrBaseRegion(CHR1, 1095, 1499);
        ChrBaseRegion hint = new ChrBaseRegion(CHR1, 1093, 1497);
        Candidate cand = new Candidate(
                CHR1, true, 151, 1001, "94M57S",
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));

        Result result = defaultResolver(annotated(annotatedJunc)).resolve(cand);

        assertTrue(result.merged());
        assertEquals(1095, result.introducedIntrons().get(0).start());
    }

    @Test
    public void testMateHintPinsIntronEndWhenPrimaryIsDownstream()
    {
        // The primary carries the leading clip, so the hint pins the intron end rather than its start. Same geometry as
        // the upstream case with the roles swapped, so the merged alignment is identical.
        Candidate withoutHint = candidate(
                CHR1, true, 151, 1498, "92S59M",
                supp(0, CHR1, true, 1001, "94M57S", 60));
        Result resNoHint = defaultResolver(annotated()).resolve(withoutHint);
        assertTrue(resNoHint.merged());
        assertEquals(1498, resNoHint.introducedIntrons().get(0).end());

        ChrBaseRegion hint = new ChrBaseRegion(CHR1, 1093, 1497);
        Candidate withHint = new Candidate(
                CHR1, true, 151, 1498, "92S59M",
                Arrays.asList(supp(0, CHR1, true, 1001, "94M57S", 60)),
                null, Arrays.asList(hint));
        Result resHinted = defaultResolver(annotated()).resolve(withHint);
        assertTrue(resHinted.merged());
        assertEquals(1497, resHinted.introducedIntrons().get(0).end());
        assertEquals("92M405N59M", resHinted.mergedCigar());
        assertEquals(1001, resHinted.mergedStart());
    }

    @Test
    public void testMateHintOutsideOverlapWindowIgnored()
    {
        // Falls back to the midpoint of the ambiguous range.
        ChrBaseRegion hint = new ChrBaseRegion(CHR1, 50000, 50100);
        Candidate cand = new Candidate(
                CHR1, true, 151, 1001, "94M57S",
                Arrays.asList(supp(0, CHR1, true, 1498, "92S59M", 60)),
                null, Arrays.asList(hint));

        Result result = defaultResolver(annotated()).resolve(cand);

        assertTrue(result.merged());
        assertEquals(1094, result.introducedIntrons().get(0).start());
    }

    @Test
    public void testRejectReasonIsNullOnSuccess()
    {
        ChrBaseRegion intron = new ChrBaseRegion(CHR1, 1095, 1499);
        Candidate cand = candidate(
                CHR1, true, READ_LEN, 1001, "94M57S",
                supp(0, CHR1, true, 1500, "94S57M", 60));

        Result result = defaultResolver(annotated(intron)).resolve(cand);

        assertTrue(result.merged());
        assertNull(result.rejectReason());
    }

    @Test
    public void testClassifiesMotifTiers()
    {
        // Strand is unknown at scan time, so GT-AG and its reverse complement CT-AC both score canonical.
        assertEquals(Tier.CANONICAL, SupplementaryResolver.motifTier(bases("GT"), bases("AG")));
        assertEquals(Tier.CANONICAL, SupplementaryResolver.motifTier(bases("CT"), bases("AC")));

        assertEquals(Tier.SEMI_CANONICAL, SupplementaryResolver.motifTier(bases("GC"), bases("AG")));
        assertEquals(Tier.SEMI_CANONICAL, SupplementaryResolver.motifTier(bases("CT"), bases("GC")));
        assertEquals(Tier.SEMI_CANONICAL, SupplementaryResolver.motifTier(bases("AT"), bases("AC")));
        assertEquals(Tier.SEMI_CANONICAL, SupplementaryResolver.motifTier(bases("GT"), bases("AT")));

        assertEquals(Tier.CANONICAL, SupplementaryResolver.motifTier(bases("gt"), bases("ag")));
        assertEquals(Tier.CANONICAL, SupplementaryResolver.motifTier(bases("ct"), bases("ac")));

        // Both flanks must match: a canonical donor with a non-motif acceptor is NONE.
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(bases("AA"), bases("GG")));
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(bases("NN"), bases("NN")));
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(bases("GT"), bases("CC")));
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(bases("CC"), bases("AG")));
    }

    @Test
    public void testRejectsNullOrWrongLengthInputs()
    {
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(null, bases("AG")));
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(bases("GT"), null));
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(bases("G"), bases("AG")));
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(bases("GT"), bases("A")));
        assertEquals(Tier.NONE, SupplementaryResolver.motifTier(bases("GTC"), bases("AG")));
    }

}
