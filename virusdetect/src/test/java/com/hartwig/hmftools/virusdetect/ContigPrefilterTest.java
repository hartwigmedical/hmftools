package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import org.junit.Test;

public class ContigPrefilterTest
{
    private static final int LENGTH = 1000;

    // The vote-density floor at this read length is 1.0 * 0.10 * 1000 / 150 ~= 0.67 votes.
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final OncologyGroup GROUP = new OncologyGroup("Group A");
    private static final ViralContig CONTIG = new ViralContig("v1", LENGTH, "Virus v1", GROUP);
    private static final ViralContig SIBLING = new ViralContig("v2", LENGTH, "Virus v2", GROUP);

    // The sibling reaches the coverage floor, establishing the group, so this contig is kept down to the relaxed floor.
    @Test
    public void testRelaxedFloorKeepsStraddlingSibling()
    {
        assertEquals(ContigFilterStatus.CANDIDATE, statusWithPresentGroup(95));
        assertEquals(ContigFilterStatus.LOW_COVERAGE, statusWithPresentGroup(89));
    }

    // No contig reached the coverage floor, so the relaxed floor never applies.
    @Test
    public void testAbsentGroupRejectsEveryContig()
    {
        Map<ViralContig, ContigFilterStatus> statuses = statuses(Map.of(CONTIG, 95, SIBLING, 95), Map.of(CONTIG, 100.0, SIBLING, 100.0));

        assertEquals(ContigFilterStatus.LOW_COVERAGE, statuses.get(CONTIG));
        assertEquals(ContigFilterStatus.LOW_COVERAGE, statuses.get(SIBLING));
    }

    @Test
    public void testVoteDensityFloorDropsCoveredContig()
    {
        assertEquals(ContigFilterStatus.CANDIDATE, status(500, 0.67));
        assertEquals(ContigFilterStatus.LOW_VOTE_DENSITY, status(500, 0.66));
    }

    // Failing both gates is reported against coverage.
    @Test
    public void testCoverageReasonTakesPrecedence()
    {
        assertEquals(ContigFilterStatus.LOW_COVERAGE, statusWithPresentGroup(50, 0.1));
    }

    private static ContigFilterStatus status(int coveredBases, double readVotes)
    {
        return statuses(Map.of(CONTIG, coveredBases), Map.of(CONTIG, readVotes)).get(CONTIG);
    }

    // The contig under test alongside a sibling with ample coverage and votes, so the group is present.
    private static ContigFilterStatus statusWithPresentGroup(int coveredBases)
    {
        return statusWithPresentGroup(coveredBases, 100);
    }

    private static ContigFilterStatus statusWithPresentGroup(int coveredBases, double readVotes)
    {
        return statuses(
                Map.of(CONTIG, coveredBases, SIBLING, LENGTH), Map.of(CONTIG, readVotes, SIBLING, 100.0)).get(CONTIG);
    }

    private static Map<ViralContig, ContigFilterStatus> statuses(
            Map<ViralContig, Integer> coveredBases, Map<ViralContig, Double> readVotes)
    {
        return ContigPrefilter.statuses(coveredBases, readVotes, MEAN_READ_LENGTH);
    }
}
