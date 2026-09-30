package com.hartwig.hmftools.viridian.detection.align;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.viridian.detection.common.ReadId;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.junit.Test;

public class ViralReadAlignmentsTest
{
    private static final int LENGTH = 100;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final ViralContig V1 = new ViralContig("v1", LENGTH, "Virus v1", GROUP_A);
    private static final ViralContig V2 = new ViralContig("v2", LENGTH, "Virus v2", GROUP_A);

    @Test
    public void testIndexesHoldEveryAlignmentBothWays()
    {
        ViralReadAlignments alignments = new ViralReadAlignments(List.of(
                alignment("r1", V1), alignment("r1", V2), alignment("r2", V1)));

        assertEquals(2, alignments.readCount());

        assertEquals(Set.of(ReadId.parse("r1"), ReadId.parse("r2")), alignments.byContig().get(V1).keySet());
        assertEquals(Set.of(ReadId.parse("r1")), alignments.byContig().get(V2).keySet());

        assertEquals(Set.of(V1, V2), alignments.byRead().get(ReadId.parse("r1")).keySet());
        assertEquals(Set.of(V1), alignments.byRead().get(ReadId.parse("r2")).keySet());
    }

    // Counting one read twice on a contig would inflate that contig's depth and read count, so the store refuses it
    // rather than leaving each caller to remember.
    @Test
    public void testRepeatAlignmentOfReadOnOneContigRejected()
    {
        List<ViralReadAlignment> repeated = List.of(alignment("r1", V1), alignment("r1", V1));

        assertThrows(IllegalArgumentException.class, () -> new ViralReadAlignments(repeated));
    }

    private static ViralReadAlignment alignment(String readName, ViralContig contig)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), contig, 1, LENGTH, LENGTH + "M", 0, 0, 100, 0, List.of(new AlignedInterval(1, LENGTH)));
    }
}
