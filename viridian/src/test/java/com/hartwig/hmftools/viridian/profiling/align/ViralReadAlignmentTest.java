package com.hartwig.hmftools.viridian.profiling.align;

import static com.hartwig.hmftools.viridian.ViridianTestUtils.readAlignment;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_GENOME_ORIGIN_CLIP_TOLERANCE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.Virus;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusGroup;

import org.junit.Test;

import htsjdk.samtools.TextCigarCodec;

public class ViralReadAlignmentTest
{
    private static final VirusGenome GENOME = new VirusGenome("v1", 200, new Virus("Virus 1", new VirusGroup("Group 1")));

    @Test
    public void testMatchedIntervals()
    {
        // Reference 30-49 aligned with a mismatch at 34, deletion at 50-51, reference 52-61 aligned with a mismatch at 55.

        ViralReadAlignment alignment = new ViralReadAlignment(
                ReadId.parse("r1"), GENOME, 30, TextCigarCodec.decode("3S20M2D10M4S"), 60, 11, new int[] { 34, 55 });

        List<MatchedInterval> expected = List.of(
                new MatchedInterval(30, 4), new MatchedInterval(35, 15), new MatchedInterval(52, 3), new MatchedInterval(56, 6));
        assertEquals(expected, alignment.matchedIntervals());
    }

    @Test
    public void testMatchedIntervalsMismatchOutsideAlignment()
    {
        ViralReadAlignment alignment = new ViralReadAlignment(
                ReadId.parse("r1"), GENOME, 30, TextCigarCodec.decode("20M2D10M"), 60, 3, new int[] { 62 });

        assertThrows(IllegalArgumentException.class, alignment::matchedIntervals);
    }

    @Test
    public void testConstructorPastGenomeEnd()
    {
        assertThrows(IllegalArgumentException.class, () -> readAlignment("r1", GENOME, 180, "30M", 60, 0));
    }

    @Test
    public void testClipsOverGenomeEnd()
    {
        // A clip projecting past a virus genome end straddles the circular origin. Up to the tolerance is allowed.
        // The 200 base genome is aligned at 11-60 with a left clip, and at 151-200 with a right clip.

        int tolerance = VIRUS_GENOME_ORIGIN_CLIP_TOLERANCE;
        assertFalse(leftClipped(10 + tolerance).clipsOverGenomeEnd());
        assertTrue(leftClipped(11 + tolerance).clipsOverGenomeEnd());
        assertFalse(rightClipped(tolerance).clipsOverGenomeEnd());
        assertTrue(rightClipped(tolerance + 1).clipsOverGenomeEnd());
    }

    @Test
    public void testBestFirstTotalOrderDivergence()
    {
        // Soft clips cost nothing in the aligner score but are bases the virus genome fails to explain, so divergence
        // and score can disagree.

        ViralReadAlignment wholeRead = readAlignment("r1", GENOME, 10, "50M", 20, 6);
        ViralReadAlignment clipped = readAlignment("r1", GENOME, 30, "20S30M", 25, 20);

        assertEquals(wholeRead, best(clipped, wholeRead));
    }

    @Test
    public void testBestFirstTotalOrderDivergenceTie()
    {
        // Equal divergence can still arise from different gap structures.

        ViralReadAlignment weaker = readAlignment("r1", GENOME, 10, "50M", 60, 8);
        ViralReadAlignment stronger = readAlignment("r1", GENOME, 10, "50M", 75, 8);

        assertEquals(stronger, best(weaker, stronger));
    }

    private static ViralReadAlignment leftClipped(int clipLength)
    {
        return readAlignment("r1", GENOME, 11, clipLength + "S50M", 50, clipLength);
    }

    private static ViralReadAlignment rightClipped(int clipLength)
    {
        return readAlignment("r1", GENOME, 151, "50M" + clipLength + "S", 50, clipLength);
    }

    private static ViralReadAlignment best(ViralReadAlignment... alignments)
    {
        return Stream.of(alignments).min(ViralReadAlignment.BEST_FIRST_TOTAL_ORDER).orElseThrow();
    }
}
