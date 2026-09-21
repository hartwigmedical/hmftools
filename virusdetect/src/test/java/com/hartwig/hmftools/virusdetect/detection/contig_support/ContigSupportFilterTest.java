package com.hartwig.hmftools.virusdetect.detection.contig_support;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.virusdetect.reference.OncologyGroup;
import com.hartwig.hmftools.virusdetect.reference.ViralContig;

import org.junit.Test;

public class ContigSupportFilterTest
{
    private static final int LENGTH = 1000;

    // The vote-density floor at this read length is 1.0 * 0.10 * 1000 / 150 ~= 0.67 votes.
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final OncologyGroup GROUP = new OncologyGroup("Group A");
    private static final ViralContig CONTIG = new ViralContig("v1", LENGTH, "Virus v1", GROUP);
    private static final ViralContig SIBLING = new ViralContig("v2", LENGTH, "Virus v2", GROUP);

    @Test
    public void testRelaxedFloorKeepsStraddlingSibling()
    {
        // The sibling establishes the group, so this contig is kept down to the relaxed floor of 90 bases.
        assertEquals(ContigFilterStatus.CANDIDATE, statusWithPresentGroup(95));
        assertEquals(ContigFilterStatus.LOW_COVERAGE, statusWithPresentGroup(89));
    }

    @Test
    public void testAbsentGroupRejectsEveryContig()
    {
        // No contig reaches the coverage floor, so the relaxed floor never applies.
        Map<ViralContig, ContigFilterStatus> statuses = statuses(Map.of(CONTIG, 95, SIBLING, 95), Map.of(CONTIG, 100.0, SIBLING, 100.0));

        assertEquals(ContigFilterStatus.LOW_COVERAGE, statuses.get(CONTIG));
        assertEquals(ContigFilterStatus.LOW_COVERAGE, statuses.get(SIBLING));
    }

    @Test
    public void testGroupPresenceAtCoverageFloor()
    {
        // A lone contig establishes its group only by reaching the coverage floor, here 0.10 * 1000 = 100 bases.
        assertEquals(ContigFilterStatus.CANDIDATE, status(100, 100.0));
        assertEquals(ContigFilterStatus.LOW_COVERAGE, status(99, 100.0));
    }

    @Test
    public void testVoteDensityFloorDropsCoveredContig()
    {
        assertEquals(ContigFilterStatus.CANDIDATE, status(500, 0.67));
        assertEquals(ContigFilterStatus.LOW_VOTE_DENSITY, status(500, 0.66));
    }

    @Test
    public void testCoverageReasonTakesPrecedence()
    {
        // Both filters fail here.
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
        return ContigSupportFilter.statuses(coveredBases, readVotes, MEAN_READ_LENGTH);
    }
}
