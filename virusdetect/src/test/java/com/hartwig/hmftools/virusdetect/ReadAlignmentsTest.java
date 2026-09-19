package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

public class ReadAlignmentsTest
{
    private static final ViralContig V1 = new ViralContig("v1", 1000, "Virus v1", new OncologyGroup("Group A"));

    // Soft clips cost nothing in the aligner score but are bases the contig fails to explain, so the two can disagree.
    // The alignment explaining more of the read wins even when the other scores higher.
    @Test
    public void testLowestDivergenceBeatsHighestScore()
    {
        ViralAlignment wholeRead = alignment(100, 6, 70);   // 100M with mismatches
        ViralAlignment clipped = alignment(100, 20, 80);    // 80M20S, cleaner but shorter

        assertEquals(wholeRead, best(List.of(clipped, wholeRead)));
    }

    // Equal divergence can still arise from different gap structures, so the stronger alignment breaks the tie.
    @Test
    public void testHighestScoreBreaksDivergenceTie()
    {
        ViralAlignment weaker = alignment(100, 8, 60);
        ViralAlignment stronger = alignment(200, 8, 75);

        assertEquals(stronger, best(List.of(weaker, stronger)));
    }

    // With nothing else to separate them, the leftmost alignment keeps the choice deterministic.
    @Test
    public void testAlignmentStartBreaksRemainingTie()
    {
        ViralAlignment later = alignment(500, 4, 90);
        ViralAlignment earlier = alignment(100, 4, 90);

        assertEquals(earlier, best(List.of(later, earlier)));
    }

    // Repeat alignments to one contig collapse to a single hit, counted so the repeat stays visible.
    @Test
    public void testRepeatAlignmentsCollapseToOneHit()
    {
        ReadAlignments read = ReadAlignments.from("r1", List.of(alignment(100, 4, 90), alignment(500, 9, 80)));

        assertEquals(1, read.hits().size());
        assertEquals(2, read.hits().get(V1).alignmentCount());
        assertEquals(4, read.hits().get(V1).divergence());
    }

    private static ViralAlignment best(List<ViralAlignment> alignments)
    {
        return ReadAlignments.from("r1", alignments).hits().get(V1).best();
    }

    private static ViralAlignment alignment(int start, int divergence, int alignerScore)
    {
        return new ViralAlignment(
                "r1", V1, start, start + 99, 0, 0, alignerScore, divergence,
                List.of(new AlignedInterval(start, 100)));
    }
}
