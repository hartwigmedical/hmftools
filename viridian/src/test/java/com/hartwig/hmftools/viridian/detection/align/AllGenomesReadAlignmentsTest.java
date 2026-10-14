package com.hartwig.hmftools.viridian.detection.align;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.common.ReadId;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class AllGenomesReadAlignmentsTest
{
    private static final int LENGTH = 100;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_H = new OncologyGroup("Group H");

    private static final VirusGenome V1 = new VirusGenome("v1", LENGTH, "Virus v1", GROUP_A);
    private static final VirusGenome V2 = new VirusGenome("v2", LENGTH, "Virus v2", GROUP_A);
    private static final VirusGenome H1 = new VirusGenome("h1", LENGTH, "Virus h1", GROUP_H);

    @Test
    public void testFromOriginStraddlersExcludedAndCounted()
    {
        // An alignment clipping over a contig end straddles the circular origin.
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                List.of(alignment("r1", V1), straddler("r2", V1), straddler("r3", V1)), 150.0);

        assertEquals(1, allAlignments.store().readCount());
        assertEquals(Map.of(V1, 2), allAlignments.metrics().originClippedReads());
    }

    @Test
    public void testFromReadCountedOncePerOncologyGroup()
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                List.of(alignment("r1", V1), alignment("r1", V2), alignment("r2", V1)), 150.0);

        assertEquals(Map.of(GROUP_A, 2), allAlignments.metrics().readCountsByOncologyGroup());
    }

    @Test
    public void testFromReadCountedInEveryOncologyGroupItAligns()
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(List.of(alignment("r1", V1), alignment("r1", H1)), 150.0);

        assertEquals(Map.of(GROUP_A, 1, GROUP_H, 1), allAlignments.metrics().readCountsByOncologyGroup());
    }

    @Test
    public void testFromStraddlersDoNotCountTowardsReadCounts()
    {
        // A contig carrying only straddlers leaves its group unrepresented.
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(List.of(alignment("r1", V1), straddler("r2", H1)), 150.0);

        assertEquals(Map.of(GROUP_A, 1), allAlignments.metrics().readCountsByOncologyGroup());
    }

    @Test
    public void testFromOriginClippedCountsReadsNotAlignments()
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                List.of(straddler("r1", V1), straddler("r1", V1), straddler("r2", V1)), 150.0);

        assertEquals(Map.of(V1, 2), allAlignments.metrics().originClippedReads());
    }

    // The aligner reports every placement it finds, so one read can be aligned to a contig several times. Only its
    // best alignment is kept, and how many it had is measured instead.
    @Test
    public void testFromRepeatAlignmentsReducedToBestAndCounted()
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                List.of(
                        divergentAlignment("r1", V1, 4),
                        divergentAlignment("r1", V1, 9),
                        divergentAlignment("r2", V1, 1)),
                150.0);

        Map<ReadId, ViralReadAlignment> contigAlignments = allAlignments.store().byContig().get(V1);
        assertEquals(2, contigAlignments.size());
        assertEquals(4, contigAlignments.get(ReadId.parse("r1")).divergence());
        assertEquals(Map.of(V1, List.of(2, 1)), allAlignments.metrics().alignmentCountsByContig());
    }

    // A sample with no viral reads has no read length to measure; this used to fail the whole run.
    @Test
    public void testFromNoAlignmentsHasNoMeanReadLength()
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(List.of(), null);

        assertEquals(0, allAlignments.store().readCount());
        assertEquals(new AllGenomesReadAlignmentsMetrics(null, Map.of(), Map.of(), Map.of()), allAlignments.metrics());
    }

    @Test
    public void testFromMeanReadLengthPresentExactlyWhenAligned()
    {
        assertThrows(IllegalArgumentException.class, () -> AllGenomesReadAlignments.from(List.of(), 150.0));
        assertThrows(IllegalArgumentException.class, () -> AllGenomesReadAlignments.from(List.of(alignment("r1", V1)), null));
        assertThrows(IllegalArgumentException.class, () -> AllGenomesReadAlignments.from(List.of(straddler("r1", V1)), null));
    }

    private static ViralReadAlignment divergentAlignment(String readName, VirusGenome contig, int divergence)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), contig, 1, LENGTH, LENGTH + "M", 0, 0, 100 - divergence, divergence,
                List.of(new AlignedInterval(1, LENGTH)));
    }

    private static ViralReadAlignment alignment(String readName, VirusGenome contig)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), contig, 1, LENGTH, LENGTH + "M", 0, 0, 100, 0, List.of(new AlignedInterval(1, LENGTH)));
    }

    // A right clip projecting well past the contig end.
    private static ViralReadAlignment straddler(String readName, VirusGenome contig)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), contig, 1, LENGTH, LENGTH + "M50S", 0, 50, 100, 50, List.of(new AlignedInterval(1, LENGTH)));
    }
}
