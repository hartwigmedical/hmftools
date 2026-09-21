package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

public class ReadAlignmentsTest
{
    private static final int LENGTH = 100;

    private static final ViralContig V1 = new ViralContig("v1", 1000, "Virus v1", new OncologyGroup("Group A"));

    @Test
    public void testRepeatAlignmentsCollapseToOneHit()
    {
        ReadAlignments read = ReadAlignments.from("r1", List.of(alignment(100, 4, 90), alignment(500, 9, 80)));

        assertEquals(1, read.hits().size());
        assertEquals(2, read.hits().get(V1).alignmentCount());
        assertEquals(4, read.hits().get(V1).divergence());
    }

    private static ViralAlignment alignment(int start, int divergence, int alignerScore)
    {
        return new ViralAlignment(
                "r1", V1, start, start + LENGTH - 1, 0, 0, alignerScore, divergence,
                List.of(new AlignedInterval(start, LENGTH)));
    }
}
