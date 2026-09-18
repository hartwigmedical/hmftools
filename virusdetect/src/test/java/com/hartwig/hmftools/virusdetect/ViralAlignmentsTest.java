package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;
import java.util.Map;

import org.junit.Test;

public class ViralAlignmentsTest
{
    private static final int LENGTH = 100;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_H = new OncologyGroup("Group H");

    private static final ViralContig V1 = new ViralContig("v1", LENGTH, "Virus v1", GROUP_A);
    private static final ViralContig V2 = new ViralContig("v2", LENGTH, "Virus v2", GROUP_A);
    private static final ViralContig H1 = new ViralContig("h1", LENGTH, "Virus h1", GROUP_H);

    // An alignment clipping over a contig end straddles the circular origin: excluded once here, and reported as a
    // per-contig count so the drop stays visible.
    @Test
    public void testOriginStraddlersExcludedAndCounted()
    {
        ViralAlignments alignments = ViralAlignments.from(
                List.of(alignment("r1", V1), straddler("r2", V1), straddler("r3", V1)), 150.0);

        assertEquals(1, alignments.alignments().size());
        assertEquals(Map.of(V1, 2), alignments.originClippedReads());
    }

    // A read aligning to several contigs of a group counts once for that group, unlike the per-contig read counts.
    @Test
    public void testReadCountedOncePerOncologyGroup()
    {
        ViralAlignments alignments = ViralAlignments.from(
                List.of(alignment("r1", V1), alignment("r1", V2), alignment("r2", V1)), 150.0);

        assertEquals(Map.of(GROUP_A, 2), alignments.readCountsByOncologyGroup());
    }

    // A read aligning across groups counts towards each of them.
    @Test
    public void testReadCountedInEveryOncologyGroupItAligns()
    {
        ViralAlignments alignments = ViralAlignments.from(List.of(alignment("r1", V1), alignment("r1", H1)), 150.0);

        assertEquals(Map.of(GROUP_A, 1, GROUP_H, 1), alignments.readCountsByOncologyGroup());
    }

    // A straddler contributes to no read count, so a contig carrying only straddlers leaves its group unrepresented.
    @Test
    public void testStraddlersDoNotCountTowardsReadCounts()
    {
        ViralAlignments alignments = ViralAlignments.from(List.of(alignment("r1", V1), straddler("r2", H1)), 150.0);

        assertEquals(Map.of(GROUP_A, 1), alignments.readCountsByOncologyGroup());
    }

    private static ViralAlignment alignment(String readName, ViralContig contig)
    {
        return new ViralAlignment(
                readName, contig, 1, LENGTH, 0, 0, 100, 0, List.of(new ViralAlignment.AlignedInterval(1, LENGTH)));
    }

    // A right clip projecting well past the contig end.
    private static ViralAlignment straddler(String readName, ViralContig contig)
    {
        return new ViralAlignment(
                readName, contig, 1, LENGTH, 0, 50, 100, 50, List.of(new ViralAlignment.AlignedInterval(1, LENGTH)));
    }
}
