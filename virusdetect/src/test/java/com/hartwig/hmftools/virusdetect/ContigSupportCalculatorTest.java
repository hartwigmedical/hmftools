package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.virusdetect.VirusConstants.READ_VOTE_CORRECT_BASE_PROBABILITY;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

public class ContigSupportCalculatorTest
{
    private static final double EPSILON = 1e-9;

    // Only the vote-density prefilter reads this; the stats themselves are unaffected.
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final ViralContig V1 = new ViralContig("v1", 20, "Virus 1", new OncologyGroup("Group 1"));
    private static final ViralContig V2 = new ViralContig("v2", 10, "Virus 2", new OncologyGroup("Group 2"));

    @Test
    public void testComputesPerContigDepthAndCoverage()
    {
        List<ViralAlignment> alignments = List.of(
                alignment("r1", V1, 1, 10, 10, 0),
                alignment("r1", V1, 1, 10, 5, 0),   // second, lower-scoring alignment of r1 on v1
                alignment("r2", V1, 6, 5, 8, 0),
                alignment("r3", V2, 1, 10, 9, 0));

        List<ContigSupport> stats = new ContigSupportCalculator(READ_VOTE_CORRECT_BASE_PROBABILITY).compute(
                ViralAlignments.from(alignments, MEAN_READ_LENGTH));

        assertEquals(2, stats.size());

        // Depth over v1's 20 bases: 1-5 from r1 only, 6-10 from both reads, 11-20 uncovered.
        int[] v1Depth = { 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        ContigSupport expectedV1 = new ContigSupport(
                V1, ContigFilterStatus.CANDIDATE, 2, 1, SummaryStats.from(new int[] { 2, 1 }), 0, 10,
                SummaryStats.from(v1Depth), SummaryStats.from(new int[] { 10, 8 }), 2.0);
        assertEquals(expectedV1, get(stats, V1));

        int[] v2Depth = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        ContigSupport expectedV2 = new ContigSupport(
                V2, ContigFilterStatus.CANDIDATE, 1, 0, SummaryStats.from(new int[] { 1 }), 0, 10,
                SummaryStats.from(v2Depth), SummaryStats.from(new int[] { 9 }), 1.0);
        assertEquals(expectedV2, get(stats, V2));
    }

    @Test
    public void testVotesAttributeVirusStrainSupport()
    {
        List<ViralAlignment> alignments = List.of(
                alignment("r1", V1, 1, 10, 10, 1),
                alignment("r1", V2, 1, 10, 6, 4),   // Aligns to another contig with more divergence
                alignment("r2", V1, 1, 10, 10, 0),
                alignment("r2", V2, 1, 10, 5, 5));  // Aligns to another contig with more divergence

        // Injected correct-base probability 0.5, so each extra divergent base halves a contig's weight (0.5^diff).
        List<ContigSupport> stats = new ContigSupportCalculator(0.5).compute(
                ViralAlignments.from(alignments, MEAN_READ_LENGTH));

        // r1: v1 weight 0.5^0, v2 0.5^3; r2: v1 0.5^0, v2 0.5^5. Votes per read sum to 1, so both contigs sum to 2.
        double v1Votes = get(stats, V1).readVotes();
        assertEquals(1.0 / (1 + Math.pow(0.5, 3)) + 1.0 / (1 + Math.pow(0.5, 5)), v1Votes, EPSILON);
        assertEquals(2.0 - v1Votes, get(stats, V2).readVotes(), EPSILON);
    }

    @Test
    public void testTiedReadSplitsVoteEvenly()
    {
        List<ViralAlignment> alignments = List.of(
                alignment("r", V1, 1, 10, 10, 2),
                alignment("r", V2, 1, 10, 10, 2));

        List<ContigSupport> stats = new ContigSupportCalculator(0.5).compute(
                ViralAlignments.from(alignments, MEAN_READ_LENGTH));

        assertEquals(0.5, get(stats, V1).readVotes(), EPSILON);
        assertEquals(0.5, get(stats, V2).readVotes(), EPSILON);
    }

    @Test
    public void testDropsAlignmentsClippingOverContigEnds()
    {
        List<ViralAlignment> alignments = List.of(
                clipped("r1", V1, 1, 10, 10, 0),    // left clip projects to -9: over the start
                clipped("r2", V1, 11, 20, 0, 10),   // right clip projects to 30: over the end
                alignment("r3", V1, 6, 10, 9, 0),   // no clip: kept
                clipped("r4", V1, 10, 19, 5, 0));   // clip projects to 5, within the contig: kept

        List<ContigSupport> stats = new ContigSupportCalculator(READ_VOTE_CORRECT_BASE_PROBABILITY).compute(
                ViralAlignments.from(alignments, MEAN_READ_LENGTH));

        ContigSupport v1 = get(stats, V1);
        assertEquals(2, v1.readCount());            // r3 and r4 kept
        assertEquals(2, v1.originClippedReads());   // r1 and r2 dropped
    }

    @Test
    public void testContigWithOnlyOriginClippedAlignmentsIsStillReported()
    {
        List<ViralAlignment> alignments = List.of(
                clipped("r1", V1, 1, 10, 30, 0),    // Clipped over origin; dropped
                clipped("r2", V1, 1, 10, 30, 0),    // Clipped over origin; dropped
                alignment("r3", V2, 1, 10, 9, 0));  // Not clipped; kept

        List<ContigSupport> stats = new ContigSupportCalculator(READ_VOTE_CORRECT_BASE_PROBABILITY).compute(
                ViralAlignments.from(alignments, MEAN_READ_LENGTH));

        ContigSupport expected = new ContigSupport(
                V1, ContigFilterStatus.LOW_COVERAGE, 0, 0, null, 2, 0, SummaryStats.from(new int[20]), null, 0.0);
        assertEquals(expected, get(stats, V1));
    }

    private static ContigSupport get(List<ContigSupport> stats, ViralContig contig)
    {
        return stats.stream().filter(support -> support.contig().equals(contig)).findFirst().orElseThrow();
    }

    // A clean, full-length alignment covering [start, start + length - 1] with no clips.
    private static ViralAlignment alignment(String readName, ViralContig contig, int start, int length, int alignerScore,
            int divergence)
    {
        return new ViralAlignment(
                readName, contig, start, start + length - 1, 0, 0, alignerScore, divergence,
                List.of(new AlignedInterval(start, length)));
    }

    // An alignment covering [start, end] with the given clip lengths hanging off each side.
    private static ViralAlignment clipped(String readName, ViralContig contig, int start, int end, int leftClip, int rightClip)
    {
        return new ViralAlignment(
                readName, contig, start, end, leftClip, rightClip, 10, leftClip + rightClip,
                List.of(new AlignedInterval(start, end - start + 1)));
    }
}
