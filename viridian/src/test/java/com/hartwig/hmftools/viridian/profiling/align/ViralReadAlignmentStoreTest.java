package com.hartwig.hmftools.viridian.profiling.align;

import static com.hartwig.hmftools.viridian.ViridianTestUtils.readAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.junit.Test;

public class ViralReadAlignmentStoreTest
{
    private static final VirusGroup GROUP_A = new VirusGroup("Group A");
    private static final VirusGenome V1 = new VirusGenome("v1", 100, new Virus("Virus v1", GROUP_A));
    private static final VirusGenome V2 = new VirusGenome("v2", 100, new Virus("Virus v2", GROUP_A));

    private static final ReadId R1 = ReadId.parse("r1");
    private static final ReadId R2 = ReadId.parse("r2");

    @Test
    public void testIndexesHoldEveryAlignmentBothWays()
    {
        ViralReadAlignment r1OnV1 = alignment("r1", V1);
        ViralReadAlignment r1OnV2 = alignment("r1", V2);
        ViralReadAlignment r2OnV1 = alignment("r2", V1);

        ViralReadAlignmentStore alignments = new ViralReadAlignmentStore(List.of(r1OnV1, r1OnV2, r2OnV1));

        assertEquals(Map.of(R1, Map.of(V1, r1OnV1, V2, r1OnV2), R2, Map.of(V1, r2OnV1)), alignments.byRead());
        assertEquals(Map.of(V1, Map.of(R1, r1OnV1, R2, r2OnV1), V2, Map.of(R1, r1OnV2)), alignments.byGenome());
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
        return readAlignment(readName, genome, 1, "100M", 100, 0);
    }
}
