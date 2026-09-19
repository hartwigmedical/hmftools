package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ContigPrefilterTest
{
    private static final int LENGTH = 1000;

    // The vote-density floor at this read length is 1.0 * 0.10 * 1000 / 150 ~= 0.67 votes.
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final ViralContig CONTIG = new ViralContig("v1", LENGTH, "Virus v1", new OncologyGroup("Group A"));

    @Test
    public void testCoverageEstablishesGroupPresence()
    {
        assertTrue(ContigPrefilter.establishesGroupPresence(0.10));
        assertFalse(ContigPrefilter.establishesGroupPresence(0.095));
    }

    // Once a sibling has established the group, this contig is kept down to the relaxed floor.
    @Test
    public void testRelaxedFloorKeepsStraddlingSibling()
    {
        assertEquals(ContigFilterStatus.CANDIDATE, status(0.095, 100));
        assertEquals(ContigFilterStatus.LOW_COVERAGE, status(0.089, 100));
    }

    // No contig reached the coverage floor, so the relaxed floor never applies.
    @Test
    public void testAbsentGroupRejectsEveryContig()
    {
        assertEquals(ContigFilterStatus.LOW_COVERAGE, ContigPrefilter.status(CONTIG, 0.095, 100, false, MEAN_READ_LENGTH));
    }

    @Test
    public void testVoteDensityFloorDropsCoveredContig()
    {
        assertEquals(ContigFilterStatus.CANDIDATE, status(0.5, 0.67));
        assertEquals(ContigFilterStatus.LOW_VOTE_DENSITY, status(0.5, 0.66));
    }

    // Failing both gates is reported against coverage.
    @Test
    public void testCoverageReasonTakesPrecedence()
    {
        assertEquals(ContigFilterStatus.LOW_COVERAGE, status(0.05, 0.1));
    }

    private static ContigFilterStatus status(double coverageFraction, double readVotes)
    {
        return ContigPrefilter.status(CONTIG, coverageFraction, readVotes, true, MEAN_READ_LENGTH);
    }
}
