package com.hartwig.hmftools.viridian.detection.common;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.align.AlignedInterval;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignments;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Test;

public class ContigStatsCalculatorTest
{
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final ViralContig CONTIG = new ViralContig("v1", 20, "Virus 1", new OncologyGroup("Group 1"));

    // A read spanning a deletion aligns in two blocks: the skipped positions are not evidence the contig is there.
    @Test
    public void testCalculateDeletionSpannedPositionsAreNotCovered()
    {
        ViralReadAlignment read = new ViralReadAlignment(
                "r1/1", CONTIG, 1, 15, "8M2D5M", 0, 0, 40, 2,
                List.of(new AlignedInterval(1, 8), new AlignedInterval(11, 5)));

        int[] depth = { 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 };
        ContigStats expected = new ContigStats(
                CONTIG, 1, 0, 13, SummaryStats.from(depth), SummaryStats.from(new int[] { 40 }));

        assertEquals(Map.of(CONTIG, expected), calculate(List.of(read)));
    }

    // Every alignment to this contig was dropped for straddling the origin, so only that loss is left to report.
    @Test
    public void testCalculateContigWithOnlyOriginClippedAlignments()
    {
        List<ViralReadAlignment> alignments = List.of(
                originClipped("r1/1"),
                originClipped("r2/1"));

        ContigStats expected = new ContigStats(CONTIG, 0, 2, 0, SummaryStats.from(new int[20]), null);

        assertEquals(Map.of(CONTIG, expected), calculate(alignments));
    }

    private static Map<ViralContig, ContigStats> calculate(List<ViralReadAlignment> alignments)
    {
        return ContigStatsCalculator.calculate(ViralReadAlignments.from(alignments, MEAN_READ_LENGTH));
    }

    // Clipped bases project past the contig start, so this alignment crosses the circular genome's origin.
    private static ViralReadAlignment originClipped(String readName)
    {
        return new ViralReadAlignment(
                readName, CONTIG, 1, 10, "30S10M", 30, 0, 20, 30, List.of(new AlignedInterval(1, 10)));
    }
}
