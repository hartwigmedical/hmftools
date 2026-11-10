package com.hartwig.hmftools.viridian.profiling.support;

import static com.hartwig.hmftools.viridian.ViridianTestUtils.readAlignment;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.READ_VOTE_CORRECT_BASE_PROBABILITY;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.viridian.profiling.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.junit.Test;

public class GenomeSupportCalculatorTest
{
    private static final double EPSILON = 1e-9;

    private static final VirusGenome V1 = new VirusGenome("v1", 20, new Virus("Virus 1", new VirusGroup("Group 1")));
    private static final VirusGenome V2 = new VirusGenome("v2", 10, new Virus("Virus 2", new VirusGroup("Group 2")));

    @Test
    public void testCalculatePerGenomeDepthAndCoverage()
    {
        List<ViralReadAlignment> alignments = List.of(
                readAlignment("r1", V1, 1, "10M", 10, 0),
                readAlignment("r1", V1, 1, "10M", 5, 0),   // second, lower-scoring alignment of r1 on v1
                readAlignment("r2", V1, 6, "5M", 8, 0),
                readAlignment("r3", V2, 1, "10M", 9, 0));

        List<GenomeSupport> supports = calculate(alignments, READ_VOTE_CORRECT_BASE_PROBABILITY);

        assertEquals(2, supports.size());

        // Depth over v1's 20 bases: 1-5 from r1 only, 6-10 from both reads, 11-20 uncovered.
        int[] v1Depth = { 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        GenomeSupport expectedV1 = new GenomeSupport(
                new GenomeStats(V1, 2, 0, 10, SummaryStats.from(v1Depth), SummaryStats.from(new int[] { 10, 8 })),
                1, SummaryStats.from(new int[] { 2, 1 }), 2.0, 2);
        assertEquals(expectedV1, supportOf(supports, V1));

        int[] v2Depth = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        GenomeSupport expectedV2 = new GenomeSupport(
                new GenomeStats(V2, 1, 0, 10, SummaryStats.from(v2Depth), SummaryStats.from(new int[] { 9 })),
                0, SummaryStats.from(new int[] { 1 }), 1.0, 1);
        assertEquals(expectedV2, supportOf(supports, V2));
    }

    @Test
    public void testCalculateVotesAttributeVirusStrainSupport()
    {
        List<ViralReadAlignment> alignments = List.of(
                readAlignment("r1", V1, 1, "10M", 10, 1),
                readAlignment("r1", V2, 1, "10M", 6, 4),   // Aligns to another virus genome with more divergence
                readAlignment("r2", V1, 1, "10M", 10, 0),
                readAlignment("r2", V2, 1, "10M", 5, 5));  // Aligns to another virus genome with more divergence

        // Injected correct-base probability 0.5, so each extra divergent base halves a virus genome's weight
        // (0.5^diff).
        List<GenomeSupport> supports = calculate(alignments, 0.5);

        // r1: v1 weight 0.5^0, v2 0.5^3; r2: v1 0.5^0, v2 0.5^5. Votes per read sum to 1, so both virus genomes sum
        // to 2.
        double v1Votes = supportOf(supports, V1).readVotes();
        assertEquals(1.0 / (1 + Math.pow(0.5, 3)) + 1.0 / (1 + Math.pow(0.5, 5)), v1Votes, EPSILON);
        assertEquals(2.0 - v1Votes, supportOf(supports, V2).readVotes(), EPSILON);
    }

    @Test
    public void testCalculateTiedReadSplitsVoteEvenly()
    {
        List<ViralReadAlignment> alignments = List.of(
                readAlignment("r", V1, 1, "10M", 10, 2),
                readAlignment("r", V2, 1, "10M", 10, 2));

        List<GenomeSupport> supports = calculate(alignments, 0.5);

        assertEquals(0.5, supportOf(supports, V1).readVotes(), EPSILON);
        assertEquals(0.5, supportOf(supports, V2).readVotes(), EPSILON);
    }

    @Test
    public void testCalculateTopMatchReadsCountsTiedReadForEachGenome()
    {
        VirusGenome v1b = new VirusGenome("v1b", 20, new Virus("Virus 1B", V1.virusGroup()));
        List<ViralReadAlignment> alignments = List.of(
                readAlignment("r1", V1, 1, "10M", 10, 0),
                readAlignment("r1", v1b, 1, "10M", 10, 0),    // Tied, so counts for both
                readAlignment("r2", V1, 1, "10M", 10, 0),
                readAlignment("r2", v1b, 1, "10M", 8, 2),     // v1 fits better
                readAlignment("r3", v1b, 1, "10M", 10, 0),    // Only on v1b
                readAlignment("r4", V1, 1, "10M", 9, 1),
                readAlignment("r4", V2, 1, "10M", 9, 1),      // Tied across virus groups, so counts for both
                readAlignment("r5", V1, 1, "10M", 7, 3),
                readAlignment("r5", V2, 1, "10M", 10, 0));    // v2 fits better

        List<GenomeSupport> supports = calculate(alignments, READ_VOTE_CORRECT_BASE_PROBABILITY);

        // v1 is the top match for r1, r2 and r4. v1b for r1 and r3. v2 for r4 and r5.
        assertEquals(3, supportOf(supports, V1).topMatchReads());
        assertEquals(2, supportOf(supports, v1b).topMatchReads());
        assertEquals(2, supportOf(supports, V2).topMatchReads());
    }

    @Test
    public void testCalculateGenomeWithOnlyOriginClippedAlignmentsIsStillReported()
    {
        List<ViralReadAlignment> alignments = List.of(
                readAlignment("r1", V1, 1, "30S10M", 10, 30),    // Clipped over origin; dropped
                readAlignment("r2", V1, 1, "30S10M", 10, 30),    // Clipped over origin; dropped
                readAlignment("r3", V2, 1, "10M", 9, 0));        // Not clipped; kept

        List<GenomeSupport> supports = calculate(alignments, READ_VOTE_CORRECT_BASE_PROBABILITY);

        GenomeSupport expected = new GenomeSupport(
                new GenomeStats(V1, 0, 2, 0, SummaryStats.from(new int[20]), null), 0, null, 0.0, 0);
        assertEquals(expected, supportOf(supports, V1));
    }

    private static List<GenomeSupport> calculate(List<ViralReadAlignment> alignments, double correctBaseProbability)
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(alignments.stream());
        return new GenomeSupportCalculator(correctBaseProbability).calculate(allAlignments);
    }

    private static GenomeSupport supportOf(List<GenomeSupport> supports, VirusGenome genome)
    {
        return supports.stream().filter(support -> support.genome().equals(genome)).findFirst().orElseThrow();
    }
}
