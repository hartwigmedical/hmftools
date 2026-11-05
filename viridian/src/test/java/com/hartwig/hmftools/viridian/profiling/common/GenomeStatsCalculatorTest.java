package com.hartwig.hmftools.viridian.profiling.common;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.junit.Test;

import htsjdk.samtools.TextCigarCodec;

public class GenomeStatsCalculatorTest
{
    private static final VirusGenome GENOME =
            new VirusGenome("v1", 20, new Virus("Virus 1", new VirusGroup("Group 1")));

    // A read spanning a deletion aligns in two blocks: the skipped positions are not evidence the virus genome is
    // there.
    @Test
    public void testCalculateDeletionSpannedPositionsAreNotCovered()
    {
        ViralReadAlignment read = new ViralReadAlignment(
                ReadId.parse("r1/1"), GENOME, 1, TextCigarCodec.decode("8M2D5M"), 40, 2, new int[0]);

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
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(alignments.stream());
        return GenomeStatsCalculator.calculate(
                allAlignments.store().byGenome(), allAlignments.metrics().originClippedReads());
    }

    // Clipped bases project past the virus genome start, so this alignment crosses the circular genome's origin.
    private static ViralReadAlignment originClipped(String readName)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), GENOME, 1, TextCigarCodec.decode("30S10M"), 20, 30, new int[0]);
    }
}
