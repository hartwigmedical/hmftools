package com.hartwig.hmftools.viridian.detection.support;

import static java.util.stream.Collectors.toMap;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.viridian.detection.common.GenomeStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class GenomeSupportFilterTest
{
    private static final int LENGTH = 1000;

    // The vote-density floor at this read length is 1.0 * 0.10 * 1000 / 150 ~= 0.67 votes.
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final OncologyGroup GROUP = new OncologyGroup("Group A");
    private static final VirusGenome CONTIG = new VirusGenome("v1", LENGTH, "Virus v1", GROUP);
    private static final VirusGenome SIBLING = new VirusGenome("v2", LENGTH, "Virus v2", GROUP);

    @Test
    public void testStatusesRelaxedFloorKeepsStraddlingSibling()
    {
        // The sibling reaches the coverage minimum, so this contig is kept down to the relaxed floor of 90 bases.
        assertEquals(GenomeFilterStatus.CANDIDATE, statusWithCoveredSibling(95));
        assertEquals(GenomeFilterStatus.LOW_COVERAGE, statusWithCoveredSibling(89));
    }

    @Test
    public void testStatusesAbsentGroupRejectsEveryContig()
    {
        // No contig reaches the coverage floor, so the relaxed floor never applies.
        Map<VirusGenome, GenomeFilterStatus> statuses = statuses(Map.of(CONTIG, 95, SIBLING, 95), Map.of(CONTIG, 100.0, SIBLING, 100.0));

        assertEquals(GenomeFilterStatus.LOW_COVERAGE, statuses.get(CONTIG));
        assertEquals(GenomeFilterStatus.LOW_COVERAGE, statuses.get(SIBLING));
    }

    @Test
    public void testStatusesGroupPresenceAtCoverageFloor()
    {
        // A lone contig establishes its group only by reaching the coverage floor, here 0.10 * 1000 = 100 bases.
        assertEquals(GenomeFilterStatus.CANDIDATE, status(100, 100.0));
        assertEquals(GenomeFilterStatus.LOW_COVERAGE, status(99, 100.0));
    }

    @Test
    public void testStatusesVoteDensityFloorDropsCoveredContig()
    {
        assertEquals(GenomeFilterStatus.CANDIDATE, status(500, 0.67));
        assertEquals(GenomeFilterStatus.LOW_VOTE_DENSITY, status(500, 0.66));
    }

    @Test
    public void testStatusesCoverageReasonTakesPrecedence()
    {
        // Both filters fail here.
        assertEquals(GenomeFilterStatus.LOW_COVERAGE, statusWithCoveredSibling(50, 0.1));
    }

    private static GenomeFilterStatus status(int coveredBases, double readVotes)
    {
        return statuses(Map.of(CONTIG, coveredBases), Map.of(CONTIG, readVotes)).get(CONTIG);
    }

    // The contig under test alongside a sibling with ample coverage and votes, so the group reaches the coverage
    // minimum.
    private static GenomeFilterStatus statusWithCoveredSibling(int coveredBases)
    {
        return statusWithCoveredSibling(coveredBases, 100);
    }

    private static GenomeFilterStatus statusWithCoveredSibling(int coveredBases, double readVotes)
    {
        return statuses(
                Map.of(CONTIG, coveredBases, SIBLING, LENGTH), Map.of(CONTIG, readVotes, SIBLING, 100.0)).get(CONTIG);
    }

    private static Map<VirusGenome, GenomeFilterStatus> statuses(
            Map<VirusGenome, Integer> coveredBases, Map<VirusGenome, Double> readVotes)
    {
        Map<VirusGenome, GenomeStats> stats = coveredBases.entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> contigStats(entry.getKey(), entry.getValue())));

        return GenomeSupportFilter.statuses(stats, readVotes, MEAN_READ_LENGTH);
    }

    // Only coverage matters to the filter; the remaining stats are filler.
    private static GenomeStats contigStats(VirusGenome contig, int coveredBases)
    {
        SummaryStats filler = SummaryStats.from(new int[] { 1 });
        return new GenomeStats(contig, 100, 0, coveredBases, filler, filler);
    }
}
