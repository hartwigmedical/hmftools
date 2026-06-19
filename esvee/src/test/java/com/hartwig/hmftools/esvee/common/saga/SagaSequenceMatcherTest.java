package com.hartwig.hmftools.esvee.common.saga;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsFromStr;
import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.esvee.MockBwaMemAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;

public class SagaSequenceMatcherTest
{
    private static final SagaSequenceMatcherConfig CONFIG = new SagaSequenceMatcherConfig(
            0.8,
            60,
            0.8,
            24,
            16,
            24,
            2
    );

    private static final int REF_ID_1 = 1;
    private static final int REF_ID_2 = 2;

    // ----- Helpers -----

    private static SagaAssembly deletionAssembly(final String id, final String sequence, int junctionOffset)
    {
        SagaVariant variant = new SagaVariant(
                id,
                new SagaBreakend(new BasePosition("chr1", 100), FORWARD),
                new SagaBreakend(new BasePosition("chr1", 200), REVERSE),
                ""
        );
        return new SagaAssembly(id, variant, List.of(junctionOffset), sequence);
    }

    private static SagaAssembly insertionAssembly(final String id, final String sequence, int junctionOffset1, int junctionOffset2)
    {
        // Ensure the assembly has enough bases to match; otherwise it will never match.
        assertTrue(junctionOffset1 > CONFIG.junctionOverlapMin());
        assertTrue(junctionOffset2 < sequence.length() - CONFIG.junctionOverlapMin());
        String insertSequence = sequence.substring(junctionOffset1, junctionOffset2);
        SagaVariant variant = new SagaVariant(
                id,
                new SagaBreakend(new BasePosition("chr1", 100), FORWARD),
                new SagaBreakend(new BasePosition("chr1", 101), REVERSE),
                insertSequence
        );
        return new SagaAssembly(id, variant, List.of(junctionOffset1, junctionOffset2), sequence);
    }

    private static String seq(int length)
    {
        // Sequence doesn't matter because we mock all the alignments.
        return "A".repeat(length);
    }

    private static BwaMemAlignment alignment(int refId, int refStart, int refEnd, final Orientation orient, final String cigar)
    {
        int flags = orient.isForward() ? 0 : SAMFlag.READ_REVERSE_STRAND.intValue();
        return new BwaMemAlignment(flags, refId, refStart, refEnd, 0, 0, 60, 0, alignScoreFromCigar(cigar), 0, toBwaCigar(cigar), "", "", 0, 0, 0);
    }

    private static int alignScoreFromCigar(final String cigarStr)
    {
        List<CigarElement> cigar = cigarElementsFromStr(cigarStr);
        assertFalse(cigar.stream().anyMatch(e -> e.getOperator() == CigarOperator.M));
        int matchScore = 1 * cigar.stream().filter(e -> e.getOperator() == CigarOperator.EQ).mapToInt(CigarElement::getLength).sum();
        int mismatchScore = -4 * cigar.stream().filter(e -> e.getOperator() == CigarOperator.X).mapToInt(CigarElement::getLength).sum();
        long gapOpenScore = -6 * cigar.stream().filter(e -> e.getOperator().isIndel()).count();
        long gapExtendScore = -1 * cigar.stream().filter(e -> e.getOperator().isIndel()).mapToInt(CigarElement::getLength).sum();
        long clipScore = -5 * cigar.stream().filter(e -> e.getOperator().isClipping()).count();
        return (int) (matchScore + mismatchScore + gapOpenScore + gapExtendScore + clipScore);
    }

    private static String toBwaCigar(final String cigar)
    {
        // We use = and X for calculating alignment score, but BWA never produces these.
        return cigar.replaceAll("=", "M").replaceAll("X", "M");
    }

    private static SagaSequenceMatcher matcher(final byte[] query, final List<BwaMemAlignment> alignments,
            final Map<Integer, SagaAssembly> assembliesByContigId)
    {
        MockBwaMemAligner aligner = new MockBwaMemAligner(Map.of(new String(query), alignments));
        return new SagaSequenceMatcher(CONFIG, aligner, assembliesByContigId);
    }

    // Constructs the expected SagaSequenceMatch from the same inputs given to the matcher, so we can compare the full record.
    private static SagaSequenceMatch expectedMatch(final BwaMemAlignment bwaAlignment, final byte[] query, final SagaAssembly assembly)
    {
        Cigar cigar = cigarFromStr(bwaAlignment.getCigar());
        return new SagaSequenceMatch(new SagaAlignment(bwaAlignment, cigar, query.length, assembly));
    }

    // ----- Tests: Accepted matches -----

    @Test
    public void testSingleMatchWith1QueryJunction()
    {
        // Single deletion assembly (1 SAGA junction). Full alignment → match returned.
        // Assembly: 100bp, junction at 50. Query: 100bp, junction at 50.
        // Alignment: 100=. Junction overlaps: left=50, right=50 → both >= JUNCTION_OVERLAP_MIN(24).
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(100), 50);
        BwaMemAlignment alignment = alignment(REF_ID_1, 0, 100, FORWARD, "100=");
        SagaSequenceMatcher sagaMatcher = matcher(query, List.of(alignment), Map.of(REF_ID_1, assembly));

        SagaSequenceMatch result = sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(50)), false, false);

        assertEquals(expectedMatch(alignment, query, assembly), result);
    }

    @Test
    public void testSingleMatchWith2QueryJunctions()
    {
        // Insertion assembly (2 SAGA junctions). Full alignment → match returned.
        // Assembly: 120bp, junctions at [40, 80]. Query: 120bp, junctions at [40, 80].
        // Alignment: 120=. All junction overlaps >= JUNCTION_OVERLAP_MIN(24).
        byte[] query = seq(120).getBytes();
        SagaAssembly assembly = insertionAssembly("ins1", seq(120), 40, 80);
        BwaMemAlignment alignment = alignment(REF_ID_1, 0, 120, FORWARD, "120=");
        SagaSequenceMatcher sagaMatcher = matcher(query, List.of(alignment), Map.of(REF_ID_1, assembly));

        SagaSequenceMatch result = sagaMatcher.matchBySequence(
                query, List.of(new SagaJunctionInfo(40), new SagaJunctionInfo(80)), false, false);

        assertEquals(expectedMatch(alignment, query, assembly), result);
    }

    @Test
    public void testShortQueryMatchingLongSagaAssembly()
    {
        // Query (60bp) is much shorter than the SAGA assembly (200bp), but aligns to the junction region → match allowed.
        // Alignment covers saga [70, 130): junction at 100 has left=30, right=30 ≥ JUNCTION_OVERLAP_MIN(24).
        // leftUnaligned=min(0, 70)=0, rightUnaligned=min(0, 70)=0 → comfortably passes align length filter.
        byte[] query = seq(60).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(200), 100);
        BwaMemAlignment alignment = alignment(REF_ID_1, 70, 130, FORWARD, "60=");
        SagaSequenceMatcher sagaMatcher = matcher(query, List.of(alignment), Map.of(REF_ID_1, assembly));

        SagaSequenceMatch result = sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(30)), false, false);

        assertEquals(expectedMatch(alignment, query, assembly), result);
    }

    @Test
    public void testLongQueryMatchingShortSagaAssembly()
    {
        // Query (200bp) is much longer than the SAGA assembly (100bp), aligned with large soft-clips on each side.
        // Alignment: 50S100=50S maps 100bp of the query (queryStart=50, queryEnd=150) to the full saga assembly.
        // leftUnaligned=min(50, 0)=0, rightUnaligned=min(50, 0)=0 → comfortably passes align length filter.
        // Query junction at 100: left=50, right=50 ≥ JUNCTION_OVERLAP_MIN(24).
        // Saga junction at 50: left=50, right=50 ≥ 24.
        byte[] query = seq(200).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(100), 50);
        BwaMemAlignment alignment = alignment(REF_ID_1, 0, 100, FORWARD, "50S100=50S");
        SagaSequenceMatcher sagaMatcher = matcher(query, List.of(alignment), Map.of(REF_ID_1, assembly));

        SagaSequenceMatch result = sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(100)), false, false);

        assertEquals(expectedMatch(alignment, query, assembly), result);
    }

    @Test
    public void testReverseStrandAllowed()
    {
        // Reverse strand alignment with allowReverseStrand=true → match is accepted.
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(100), 50);
        BwaMemAlignment alignment = alignment(REF_ID_1, 0, 100, REVERSE, "100=");
        SagaSequenceMatcher sagaMatcher = matcher(query, List.of(alignment), Map.of(REF_ID_1, assembly));

        SagaSequenceMatch result = sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(50)), false, true);

        assertEquals(expectedMatch(alignment, query, assembly), result);
    }

    @Test
    public void testLowerJunctionOverlap()
    {
        // Junction overlap of 20 is below JUNCTION_OVERLAP_MIN(24) but above JUNCTION_OVERLAP_MIN_LOWER(16).
        // With lowerJunctionOverlap=false: filtered. With lowerJunctionOverlap=true: accepted.
        // Query=100bp, junction at 20, alignment "100=": left=20, right=80.
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(100), 20);
        BwaMemAlignment alignment = alignment(REF_ID_1, 0, 100, FORWARD, "100=");
        SagaSequenceMatcher sagaMatcher = matcher(query, List.of(alignment), Map.of(REF_ID_1, assembly));

        List<SagaJunctionInfo> queryJunctions = List.of(new SagaJunctionInfo(20));
        assertNull("Junction overlap of 20 should be filtered with standard threshold (24)",
                sagaMatcher.matchBySequence(query, queryJunctions, false, false));
        assertEquals("Junction overlap of 20 should be accepted with lower threshold (16)",
                expectedMatch(alignment, query, assembly),
                sagaMatcher.matchBySequence(query, queryJunctions, true, false));
    }

    // ----- Tests: Multiple matches / tie-breaking -----

    @Test
    public void testMultipleMatchesTieBreakHigherScore()
    {
        // Two assemblies both pass all filters. The one with the higher alignment score should win.
        // REF_ID_1 uses "90=10S" (score=85: 90 matches - 5 clip penalty), REF_ID_2 uses "100=" (score=100).
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly1 = deletionAssembly("low_score", seq(100), 50);
        SagaAssembly assembly2 = deletionAssembly("high_score", seq(100), 50);
        BwaMemAlignment alignment1 = alignment(REF_ID_1, 0, 90, FORWARD, "90=10S");
        BwaMemAlignment alignment2 = alignment(REF_ID_2, 0, 100, FORWARD, "100=");
        SagaSequenceMatcher sagaMatcher = matcher(query, List.of(alignment1, alignment2), Map.of(REF_ID_1, assembly1, REF_ID_2, assembly2));

        SagaSequenceMatch result = sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(50)), false, false);

        assertEquals(expectedMatch(alignment2, query, assembly2), result);
    }

    @Test
    public void testMultipleMatchesTieBreakCloserAssemblyLength()
    {
        // Two assemblies, same alignment score (both "90=10S"). Tie-break: prefer the assembly whose
        // length is closer to the query length (100).
        // REF_ID_1: length=90 (distance=10). REF_ID_2: length=105 (distance=5) → REF_ID_2 wins.
        // sagaEnd=90 covers both junctions (45 and 52).
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly1 = deletionAssembly("farther_len", seq(90), 45);
        SagaAssembly assembly2 = deletionAssembly("closer_len", seq(105), 52);
        BwaMemAlignment alignment1 = alignment(REF_ID_1, 0, 90, FORWARD, "90=10S");
        BwaMemAlignment alignment2 = alignment(REF_ID_2, 0, 90, FORWARD, "90=10S");
        SagaSequenceMatcher sagaMatcher = matcher(query, List.of(alignment1, alignment2), Map.of(REF_ID_1, assembly1, REF_ID_2, assembly2));

        SagaSequenceMatch result = sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(45)), false, false);

        // assembly2 (length 105, distance 5) is closer to query length 100 than assembly1 (length 90, distance 10)
        assertEquals(expectedMatch(alignment2, query, assembly2), result);
    }

    // ----- Tests: Filtered matches -----

    @Test
    public void testFilterReverseStrand()
    {
        // Reverse strand alignment when allowReverseStrand=false → reverse_strand filter, no match.
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(100), 50);
        SagaSequenceMatcher sagaMatcher = matcher(
                query, List.of(alignment(REF_ID_1, 0, 100, REVERSE, "100=")), Map.of(REF_ID_1, assembly));

        assertNull("Reverse strand alignment should be filtered when allowReverseStrand=false",
                sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(50)), false, false));
    }

    @Test
    public void testFilterAlignScore()
    {
        // Alignment with many mismatches → low score → align_score filter, no match.
        // "50=50X": score = 50 + (-4*50) = -150. minAlignScore = ceil(100 * 0.8) = 80. -150 < 80 → fails.
        // Full alignment (no clips) and junction at 50 clearly pass all other filters.
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(100), 50);
        SagaSequenceMatcher sagaMatcher = matcher(
                query, List.of(alignment(REF_ID_1, 0, 100, FORWARD, "50=50X")), Map.of(REF_ID_1, assembly));

        assertNull("Alignment score below ratio threshold should be filtered",
                sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(50)), false, false));
    }

    @Test
    public void testFilterAlignLength()
    {
        // Only 40% of the sequence aligns → aligned_length filter, no match.
        // Query=200bp, saga=200bp, cigar="80=120S".
        // rightUnaligned=min(120, 120)=120. maxUnaligned=ceil(200*0.2)=40. 120 > 40 → fails.
        // Score (80-5=75) comfortably passes threshold (ceil(80*0.8)=64).
        byte[] query = seq(200).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(200), 40);
        SagaSequenceMatcher sagaMatcher = matcher(
                query, List.of(alignment(REF_ID_1, 0, 80, FORWARD, "80=120S")), Map.of(REF_ID_1, assembly));

        assertNull("Partial alignment covering only 40% of sequences should be filtered",
                sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(40)), false, false));
    }

    @Test
    public void testFilterQueryJunctionOverlap()
    {
        // Query junction at 90, alignment "100=" covers query [0, 100): right overlap = 100-90 = 10 < JUNCTION_OVERLAP_MIN(24).
        // The saga assembly is 200bp so rightUnaligned=min(0, 100)=0, making the align length filter clearly pass.
        // Saga junction at 50 is fully covered: left=50, right=50 ≥ 24.
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(200), 50);
        SagaSequenceMatcher sagaMatcher = matcher(
                query, List.of(alignment(REF_ID_1, 0, 100, FORWARD, "100=")), Map.of(REF_ID_1, assembly));

        assertNull("Query junction outside the aligned region should be filtered",
                sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(90)), false, false));
    }

    @Test
    public void testFilterSagaJunctionOverlap()
    {
        // Saga assembly junction at 30, but alignment only covers saga [60, 160) (sagaStart=60).
        // Left overlap: 30-60 = -30 < JUNCTION_OVERLAP_MIN(24) → saga_junction_overlap filter.
        // Query junction at 50 is fully covered: left=50, right=50 ≥ 24. Align length clearly passes.
        byte[] query = seq(100).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(160), 30);
        SagaSequenceMatcher sagaMatcher = matcher(
                query, List.of(alignment(REF_ID_1, 60, 160, FORWARD, "100=")), Map.of(REF_ID_1, assembly));

        assertNull("SAGA junction not covered by the alignment should be filtered",
                sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(50)), false, false));
    }

    @Test
    public void testFilterJunctionIndel()
    {
        // Large indel (10bp deletion) at the junction → junction_indel filter, no match.
        // Query=190bp, saga=200bp, cigar="100=10D90=": full alignment (leftUnaligned=rightUnaligned=0).
        // D element is at refStart=100, which is at/inside the junction → indelNearby.
        // Indel length 10 > JUNCTION_INDEL_MAX_LEN(2) → not ignored → junction_indel filter.
        // Score: 190-16=174, threshold=ceil(190*0.8)=152 → comfortably passes. Align length: 0 unaligned → clearly passes.
        byte[] query = seq(190).getBytes();
        SagaAssembly assembly = deletionAssembly("del1", seq(200), 100);
        SagaSequenceMatcher sagaMatcher = matcher(
                query, List.of(alignment(REF_ID_1, 0, 200, FORWARD, "100=10D90=")), Map.of(REF_ID_1, assembly));

        assertNull("Large indel near the junction should be filtered",
                sagaMatcher.matchBySequence(query, List.of(new SagaJunctionInfo(100)), false, false));
    }
}
