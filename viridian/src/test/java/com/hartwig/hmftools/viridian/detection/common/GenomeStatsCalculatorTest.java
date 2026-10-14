package com.hartwig.hmftools.viridian.detection.common;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.align.AlignedInterval;
import com.hartwig.hmftools.viridian.detection.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class GenomeStatsCalculatorTest
{
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final VirusGenome GENOME = new VirusGenome("v1", 20, "Virus 1", new OncologyGroup("Group 1"));

    // A read spanning a deletion aligns in two blocks: the skipped positions are not evidence the virus genome is
    // there.
    @Test
    public void testCalculateDeletionSpannedPositionsAreNotCovered()
    {
        ViralReadAlignment read = new ViralReadAlignment(
                ReadId.parse("r1/1"), GENOME, 1, 15, "8M2D5M", 0, 0, 40, 2,
                List.of(new AlignedInterval(1, 8), new AlignedInterval(11, 5)));

        int[] depth = { 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 };
        GenomeStats expected = new GenomeStats(
                GENOME, 1, 0, 13, SummaryStats.from(depth), SummaryStats.from(new int[] { 40 }));

        assertEquals(Map.of(GENOME, expected), calculate(List.of(read)));
    }

    // Every alignment to this virus genome was dropped for straddling the origin, so only that loss is left to report.
    @Test
    public void testCalculateGenomeWithOnlyOriginClippedAlignments()
    {
        List<ViralReadAlignment> alignments = List.of(
                originClipped("r1/1"),
                originClipped("r2/1"));

        GenomeStats expected = new GenomeStats(GENOME, 0, 2, 0, SummaryStats.from(new int[20]), null);

        assertEquals(Map.of(GENOME, expected), calculate(alignments));
    }

    private static Map<VirusGenome, GenomeStats> calculate(List<ViralReadAlignment> alignments)
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(alignments, MEAN_READ_LENGTH);
        return GenomeStatsCalculator.calculate(
                allAlignments.store().byGenome(), allAlignments.metrics().originClippedReads());
    }

    // Clipped bases project past the virus genome start, so this alignment crosses the circular genome's origin.
    private static ViralReadAlignment originClipped(String readName)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), GENOME, 1, 10, "30S10M", 30, 0, 20, 30, List.of(new AlignedInterval(1, 10)));
    }
}
