package com.hartwig.hmftools.viridian.detection.align;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.detection.common.ReadId;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Test;

public class AllAlignmentsTest
{
    private static final int LENGTH = 100;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_H = new OncologyGroup("Group H");

    private static final ViralContig V1 = new ViralContig("v1", LENGTH, "Virus v1", GROUP_A);
    private static final ViralContig V2 = new ViralContig("v2", LENGTH, "Virus v2", GROUP_A);
    private static final ViralContig H1 = new ViralContig("h1", LENGTH, "Virus h1", GROUP_H);

    @Test
    public void testFromOriginStraddlersExcludedAndCounted()
    {
        // An alignment clipping over a contig end straddles the circular origin.
        AllAlignments allAlignments = AllAlignments.from(
                List.of(alignment("r1", V1), straddler("r2", V1), straddler("r3", V1)), 150.0);

        assertEquals(1, allAlignments.alignments().readCount());
        assertEquals(Map.of(V1, 2), allAlignments.metrics().originClippedReads());
    }

    @Test
    public void testFromReadCountedOncePerOncologyGroup()
    {
        AllAlignments allAlignments = AllAlignments.from(
                List.of(alignment("r1", V1), alignment("r1", V2), alignment("r2", V1)), 150.0);

        assertEquals(Map.of(GROUP_A, 2), allAlignments.metrics().readCountsByOncologyGroup());
    }

    @Test
    public void testFromReadCountedInEveryOncologyGroupItAligns()
    {
        AllAlignments allAlignments = AllAlignments.from(List.of(alignment("r1", V1), alignment("r1", H1)), 150.0);

        assertEquals(Map.of(GROUP_A, 1, GROUP_H, 1), allAlignments.metrics().readCountsByOncologyGroup());
    }

    @Test
    public void testFromStraddlersDoNotCountTowardsReadCounts()
    {
        // A contig carrying only straddlers leaves its group unrepresented.
        AllAlignments allAlignments = AllAlignments.from(List.of(alignment("r1", V1), straddler("r2", H1)), 150.0);

        assertEquals(Map.of(GROUP_A, 1), allAlignments.metrics().readCountsByOncologyGroup());
    }

    @Test
    public void testFromOriginClippedCountsReadsNotAlignments()
    {
        AllAlignments allAlignments = AllAlignments.from(
                List.of(straddler("r1", V1), straddler("r1", V1), straddler("r2", V1)), 150.0);

        assertEquals(Map.of(V1, 2), allAlignments.metrics().originClippedReads());
    }

    // The aligner reports every placement it finds, so one read can be aligned to a contig several times. Only its
    // best alignment is kept, and how many it had is measured instead.
    @Test
    public void testFromRepeatAlignmentsReducedToBestAndCounted()
    {
        AllAlignments allAlignments = AllAlignments.from(
                List.of(
                        divergentAlignment("r1", V1, 4),
                        divergentAlignment("r1", V1, 9),
                        divergentAlignment("r2", V1, 1)),
                150.0);

        Map<ReadId, ViralReadAlignment> contigAlignments = allAlignments.alignments().byContig().get(V1);
        assertEquals(2, contigAlignments.size());
        assertEquals(4, contigAlignments.get(ReadId.parse("r1")).divergence());
        assertEquals(Map.of(V1, List.of(2, 1)), allAlignments.metrics().alignmentCountsByContig());
    }

    private static ViralReadAlignment divergentAlignment(String readName, ViralContig contig, int divergence)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), contig, 1, LENGTH, LENGTH + "M", 0, 0, 100 - divergence, divergence,
                List.of(new AlignedInterval(1, LENGTH)));
    }

    private static ViralReadAlignment alignment(String readName, ViralContig contig)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), contig, 1, LENGTH, LENGTH + "M", 0, 0, 100, 0, List.of(new AlignedInterval(1, LENGTH)));
    }

    // A right clip projecting well past the contig end.
    private static ViralReadAlignment straddler(String readName, ViralContig contig)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), contig, 1, LENGTH, LENGTH + "M50S", 0, 50, 100, 50, List.of(new AlignedInterval(1, LENGTH)));
    }
}
