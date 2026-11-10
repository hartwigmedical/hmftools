package com.hartwig.hmftools.viridian.profiling.align;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.junit.Test;

import htsjdk.samtools.TextCigarCodec;

public class AllGenomesReadAlignmentsTest
{
    private static final int LENGTH = 100;

    private static final VirusGroup GROUP_A = new VirusGroup("Group A");

    private static final VirusGenome V1 = new VirusGenome("v1", LENGTH, new Virus("Virus v1", GROUP_A));

    @Test
    public void testFromOriginStraddlersExcludedAndCounted()
    {
        // An alignment clipping over a virus genome end straddles the circular origin.
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                Stream.of(alignment("r1", V1), straddler("r2", V1), straddler("r3", V1)));

        assertEquals(Set.of(ReadId.parse("r1")), allAlignments.store().byRead().keySet());
        assertEquals(Map.of(V1, 2), allAlignments.originClippedReads());
    }

    @Test
    public void testFromOriginClippedCountsReadsNotAlignments()
    {
        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                Stream.of(straddler("r1", V1), straddler("r1", V1), straddler("r2", V1)));

        assertEquals(Map.of(V1, 2), allAlignments.originClippedReads());
    }

    @Test
    public void testFromRepeatAlignmentsReducedToBestAndCounted()
    {
        // The aligner reports every placement it finds, so one read can be aligned to a virus genome several times. Only
        // its best alignment is kept, wherever it comes in the input, and how many it had is measured instead.

        AllGenomesReadAlignments allAlignments = AllGenomesReadAlignments.from(
                Stream.of(
                        divergentAlignment("r1", V1, 9),
                        divergentAlignment("r1", V1, 4),
                        divergentAlignment("r1", V1, 7),
                        divergentAlignment("r2", V1, 1)));

        Map<ReadId, ViralReadAlignment> expected = Map.of(
                ReadId.parse("r1"), divergentAlignment("r1", V1, 4),
                ReadId.parse("r2"), divergentAlignment("r2", V1, 1));
        assertEquals(expected, allAlignments.store().byGenome().get(V1));
        assertEquals(Map.of(V1, List.of(3, 1)), allAlignments.alignmentCountsByGenome());
    }

    private static ViralReadAlignment divergentAlignment(String readName, VirusGenome genome, int divergence)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), genome, 1, TextCigarCodec.decode(LENGTH + "M"),
                100 - divergence, divergence, new int[0]);
    }

    private static ViralReadAlignment alignment(String readName, VirusGenome genome)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), genome, 1, TextCigarCodec.decode(LENGTH + "M"), 100, 0, new int[0]);
    }

    // A right clip projecting well past the virus genome end.
    private static ViralReadAlignment straddler(String readName, VirusGenome genome)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), genome, 1, TextCigarCodec.decode(LENGTH + "M50S"), 100, 50, new int[0]);
    }
}
