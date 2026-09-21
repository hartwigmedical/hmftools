package com.hartwig.hmftools.viridian.detection.read_align;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Test;

public class ViralReadAlignmentsTest
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
        ViralReadAlignments alignments = ViralReadAlignments.from(
                List.of(alignment("r1", V1), straddler("r2", V1), straddler("r3", V1)), 150.0);

        assertEquals(1, alignments.reads().size());
        assertEquals(Map.of(V1, 2), alignments.originClippedReads());
    }

    @Test
    public void testFromReadCountedOncePerOncologyGroup()
    {
        ViralReadAlignments alignments = ViralReadAlignments.from(
                List.of(alignment("r1", V1), alignment("r1", V2), alignment("r2", V1)), 150.0);

        assertEquals(Map.of(GROUP_A, 2), alignments.readCountsByOncologyGroup());
    }

    @Test
    public void testFromReadCountedInEveryOncologyGroupItAligns()
    {
        ViralReadAlignments alignments = ViralReadAlignments.from(List.of(alignment("r1", V1), alignment("r1", H1)), 150.0);

        assertEquals(Map.of(GROUP_A, 1, GROUP_H, 1), alignments.readCountsByOncologyGroup());
    }

    @Test
    public void testFromStraddlersDoNotCountTowardsReadCounts()
    {
        // A contig carrying only straddlers leaves its group unrepresented.
        ViralReadAlignments alignments = ViralReadAlignments.from(List.of(alignment("r1", V1), straddler("r2", H1)), 150.0);

        assertEquals(Map.of(GROUP_A, 1), alignments.readCountsByOncologyGroup());
    }

    @Test
    public void testFromOriginClippedCountsReadsNotAlignments()
    {
        ViralReadAlignments alignments = ViralReadAlignments.from(
                List.of(straddler("r1", V1), straddler("r1", V1), straddler("r2", V1)), 150.0);

        assertEquals(Map.of(V1, 2), alignments.originClippedReads());
    }

    private static ViralReadAlignment alignment(String readName, ViralContig contig)
    {
        return new ViralReadAlignment(
                readName, contig, 1, LENGTH, 0, 0, 100, 0, List.of(new AlignedInterval(1, LENGTH)));
    }

    // A right clip projecting well past the contig end.
    private static ViralReadAlignment straddler(String readName, ViralContig contig)
    {
        return new ViralReadAlignment(
                readName, contig, 1, LENGTH, 0, 50, 100, 50, List.of(new AlignedInterval(1, LENGTH)));
    }
}
