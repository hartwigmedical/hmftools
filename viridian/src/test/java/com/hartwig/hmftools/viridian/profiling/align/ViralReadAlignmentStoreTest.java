package com.hartwig.hmftools.viridian.profiling.align;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.junit.Test;

public class ViralReadAlignmentStoreTest
{
    private static final int LENGTH = 100;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final VirusGenome V1 = new VirusGenome("v1", LENGTH, "Virus v1", GROUP_A);
    private static final VirusGenome V2 = new VirusGenome("v2", LENGTH, "Virus v2", GROUP_A);

    @Test
    public void testIndexesHoldEveryAlignmentBothWays()
    {
        ViralReadAlignmentStore alignments = new ViralReadAlignmentStore(List.of(
                alignment("r1", V1), alignment("r1", V2), alignment("r2", V1)));

        assertEquals(Set.of(ReadId.parse("r1"), ReadId.parse("r2")), alignments.byRead().keySet());

        assertEquals(Set.of(ReadId.parse("r1"), ReadId.parse("r2")), alignments.byGenome().get(V1).keySet());
        assertEquals(Set.of(ReadId.parse("r1")), alignments.byGenome().get(V2).keySet());

        assertEquals(Set.of(V1, V2), alignments.byRead().get(ReadId.parse("r1")).keySet());
        assertEquals(Set.of(V1), alignments.byRead().get(ReadId.parse("r2")).keySet());
    }

    // Counting one read twice on a virus genome would inflate that genome's depth and read count, so the store refuses
    // it rather than leaving each caller to remember.
    @Test
    public void testRepeatAlignmentOfReadOnOneGenomeRejected()
    {
        List<ViralReadAlignment> repeated = List.of(alignment("r1", V1), alignment("r1", V1));

        assertThrows(IllegalArgumentException.class, () -> new ViralReadAlignmentStore(repeated));
    }

    private static ViralReadAlignment alignment(String readName, VirusGenome genome)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), genome, 1, LENGTH, LENGTH + "M", 0, 0, 100, 0, List.of(new AlignedInterval(1, LENGTH)));
    }
}
