package com.hartwig.hmftools.viridian.profiling.support;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.READ_VOTE_CORRECT_BASE_PROBABILITY;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.align.AlignedInterval;
import com.hartwig.hmftools.viridian.profiling.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStatsCalculator;
import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class GenomeSupportCalculatorTest
{
    private static final double EPSILON = 1e-9;

    // Only the vote-density prefilter reads this; the stats themselves are unaffected.
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final VirusGenome V1 = new VirusGenome("v1", 20, "Virus 1", new OncologyGroup("Group 1"));
    private static final VirusGenome V2 = new VirusGenome("v2", 10, "Virus 2", new OncologyGroup("Group 2"));

    @Test
    public void testComputePerGenomeDepthAndCoverage()
    {
        List<ViralReadAlignment> alignments = List.of(
                alignment("r1", V1, 1, 10, 10, 0),
                alignment("r1", V1, 1, 10, 5, 0),   // second, lower-scoring alignment of r1 on v1
                alignment("r2", V1, 6, 5, 8, 0),
                alignment("r3", V2, 1, 10, 9, 0));

        List<GenomeSupport> stats = compute(alignments, READ_VOTE_CORRECT_BASE_PROBABILITY);

        assertEquals(2, stats.size());

        // Depth over v1's 20 bases: 1-5 from r1 only, 6-10 from both reads, 11-20 uncovered.
        int[] v1Depth = { 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        GenomeSupport expectedV1 = new GenomeSupport(
                new GenomeStats(V1, 2, 0, 10, SummaryStats.from(v1Depth), SummaryStats.from(new int[] { 10, 8 })),
                GenomeFilterStatus.PASS, 1, SummaryStats.from(new int[] { 2, 1 }), 2.0);
        assertEquals(expectedV1, get(stats, V1));

        int[] v2Depth = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        GenomeSupport expectedV2 = new GenomeSupport(
                new GenomeStats(V2, 1, 0, 10, SummaryStats.from(v2Depth), SummaryStats.from(new int[] { 9 })),
                GenomeFilterStatus.PASS, 0, SummaryStats.from(new int[] { 1 }), 1.0);
        assertEquals(expectedV2, get(stats, V2));
    }

    @Test
    public void testComputeVotesAttributeVirusStrainSupport()
    {
        List<ViralReadAlignment> alignments = List.of(
                alignment("r1", V1, 1, 10, 10, 1),
                alignment("r1", V2, 1, 10, 6, 4),   // Aligns to another virus genome with more divergence
                alignment("r2", V1, 1, 10, 10, 0),
                alignment("r2", V2, 1, 10, 5, 5));  // Aligns to another virus genome with more divergence

        // Injected correct-base probability 0.5, so each extra divergent base halves a virus genome's weight
        // (0.5^diff).
        List<GenomeSupport> stats = compute(alignments, 0.5);

        // r1: v1 weight 0.5^0, v2 0.5^3; r2: v1 0.5^0, v2 0.5^5. Votes per read sum to 1, so both virus genomes sum
        // to 2.
        double v1Votes = get(stats, V1).readVotes();
        assertEquals(1.0 / (1 + Math.pow(0.5, 3)) + 1.0 / (1 + Math.pow(0.5, 5)), v1Votes, EPSILON);
        assertEquals(2.0 - v1Votes, get(stats, V2).readVotes(), EPSILON);
    }

    @Test
    public void testComputeTiedReadSplitsVoteEvenly()
    {
        List<ViralReadAlignment> alignments = List.of(
                alignment("r", V1, 1, 10, 10, 2),
                alignment("r", V2, 1, 10, 10, 2));

        List<GenomeSupport> stats = compute(alignments, 0.5);

        assertEquals(0.5, get(stats, V1).readVotes(), EPSILON);
        assertEquals(0.5, get(stats, V2).readVotes(), EPSILON);
    }

    @Test
    public void testComputeDropsAlignmentsClippingOverGenomeEnds()
    {
        List<ViralReadAlignment> alignments = List.of(
                clipped("r1", V1, 1, 10, 10, 0),    // left clip projects to -9: over the start
                clipped("r2", V1, 11, 20, 0, 10),   // right clip projects to 30: over the end
                alignment("r3", V1, 6, 10, 9, 0),   // no clip: kept
                clipped("r4", V1, 10, 19, 5, 0));   // clip projects to 5, within the virus genome: kept

        List<GenomeSupport> stats = compute(alignments, READ_VOTE_CORRECT_BASE_PROBABILITY);

        GenomeSupport v1 = get(stats, V1);
        assertEquals(2, v1.stats().readCount());            // r3 and r4 kept
        assertEquals(2, v1.stats().originClippedReads());   // r1 and r2 dropped
    }

    @Test
    public void testComputeGenomeWithOnlyOriginClippedAlignmentsIsStillReported()
    {
        List<ViralReadAlignment> alignments = List.of(
                clipped("r1", V1, 1, 10, 30, 0),    // Clipped over origin; dropped
                clipped("r2", V1, 1, 10, 30, 0),    // Clipped over origin; dropped
                alignment("r3", V2, 1, 10, 9, 0));  // Not clipped; kept

        List<GenomeSupport> stats = compute(alignments, READ_VOTE_CORRECT_BASE_PROBABILITY);

        GenomeSupport expected = new GenomeSupport(
                new GenomeStats(V1, 0, 2, 0, SummaryStats.from(new int[20]), null), GenomeFilterStatus.LOW_COVERAGE, 0, null, 0.0);
        assertEquals(expected, get(stats, V1));
    }

    private static List<GenomeSupport> compute(List<ViralReadAlignment> alignments, double correctBaseProbability)
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(alignments, MEAN_READ_LENGTH);
        Map<VirusGenome, GenomeStats> genomeStats = GenomeStatsCalculator.calculate(
                allAlignments.store().byGenome(), allAlignments.metrics().originClippedReads());
        return new GenomeSupportCalculator(correctBaseProbability).calculate(allAlignments, genomeStats);
    }

    private static GenomeSupport get(List<GenomeSupport> stats, VirusGenome genome)
    {
        return stats.stream().filter(support -> support.genome().equals(genome)).findFirst().orElseThrow();
    }

    // A clean, full-length alignment covering [start, start + length - 1] with no clips.
    private static ViralReadAlignment alignment(String readName, VirusGenome genome, int start, int length, int alignerScore,
            int divergence)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), genome, start, start + length - 1, length + "M", 0, 0, alignerScore, divergence,
                List.of(new AlignedInterval(start, length)));
    }

    // An alignment covering [start, end] with the given clip lengths hanging off each side.
    private static ViralReadAlignment clipped(String readName, VirusGenome genome, int start, int end, int leftClip, int rightClip)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), genome, start, end,
                leftClip + "S" + (end - start + 1) + "M" + rightClip + "S",
                leftClip, rightClip, 10, leftClip + rightClip,
                List.of(new AlignedInterval(start, end - start + 1)));
    }
}
