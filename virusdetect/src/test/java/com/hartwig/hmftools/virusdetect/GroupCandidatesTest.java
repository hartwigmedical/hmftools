package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

public class GroupCandidatesTest
{
    private static final int LENGTH = 1000;

    // The vote-density floor at this read length is 0.10 * 1000 / 150 ~= 0.67 votes.
    private static final double MEAN_READ_LENGTH = 150.0;

    // One contig clears the coverage floor, so its sibling is kept down to the relaxed floor.
    @Test
    public void testRelaxedFloorKeepsStraddlingSibling()
    {
        GroupCandidates candidates = GroupCandidates.prefilter(
                List.of(stats("v1", 0.5, 100), stats("v2", 0.095, 100)), MEAN_READ_LENGTH);

        assertEquals(List.of("v1", "v2"), names(candidates.candidates()));
        assertTrue(candidates.rejected().isEmpty());
    }

    @Test
    public void testRelaxedFloorDropsContigBelowIt()
    {
        GroupCandidates candidates = GroupCandidates.prefilter(
                List.of(stats("v1", 0.5, 100), stats("v2", 0.05, 100)), MEAN_READ_LENGTH);

        assertEquals(List.of("v1"), names(candidates.candidates()));
        assertEquals(ContigFilterStatus.LOW_COVERAGE, candidates.rejected().get(0).reason());
    }

    // No contig reaches the coverage floor, so the relaxed floor never applies and the group keeps nothing.
    @Test
    public void testAbsentGroupRejectsEveryContig()
    {
        GroupCandidates candidates = GroupCandidates.prefilter(
                List.of(stats("v1", 0.095, 100), stats("v2", 0.09, 100)), MEAN_READ_LENGTH);

        assertTrue(candidates.candidates().isEmpty());
        assertEquals(2, candidates.rejected().size());
    }

    @Test
    public void testVoteDensityFloorDropsCoveredContig()
    {
        GroupCandidates candidates = GroupCandidates.prefilter(
                List.of(stats("v1", 0.5, 100), stats("v2", 0.5, 0.5)), MEAN_READ_LENGTH);

        assertEquals(List.of("v1"), names(candidates.candidates()));
        assertEquals(ContigFilterStatus.LOW_VOTE_DENSITY, candidates.rejected().get(0).reason());
    }

    // Failing both gates is reported against coverage.
    @Test
    public void testCoverageReasonTakesPrecedence()
    {
        GroupCandidates candidates = GroupCandidates.prefilter(
                List.of(stats("v1", 0.5, 100), stats("v2", 0.05, 0.1)), MEAN_READ_LENGTH);

        assertEquals(ContigFilterStatus.LOW_COVERAGE, candidates.rejected().get(0).reason());
    }

    private static List<String> names(List<ContigSupport> stats)
    {
        return stats.stream().map(stat -> stat.contig().name()).toList();
    }

    private static ContigSupport stats(String contig, double coverage, double votes)
    {
        SummaryStats dummy = SummaryStats.from(new int[] { 1 });
        ViralContig viralContig = new ViralContig(contig, LENGTH, "Virus " + contig, new OncologyGroup("Group A"));
        return new ContigSupport(
                viralContig, 100, 0, dummy, 0, (int) Math.round(coverage * LENGTH), dummy, dummy, votes);
    }
}
