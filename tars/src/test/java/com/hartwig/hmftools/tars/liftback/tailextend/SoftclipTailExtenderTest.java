package com.hartwig.hmftools.tars.liftback.tailextend;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.repeatedBase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;
import com.hartwig.hmftools.tars.liftback.rescue.AnnotatedJunctionIndex;

import org.junit.Test;

public class SoftclipTailExtenderTest
{
    private static final String CHR1 = "chr1";

    private static TestGenome genome()
    {
        return new TestGenome().with(CHR1, 200, 'A');
    }

    private static SoftclipTailExtender extender(final TestGenome genome, final Set<ChrBaseRegion> introns)
    {
        return new SoftclipTailExtender(
                genome.asRefSource(),
                introns == null ? null : new AnnotatedJunctionIndex(introns),
                TailExtensionConfig.enabledDefaults());
    }

    private static SoftclipTailExtender extender(final TestGenome genome)
    {
        return extender(genome, null);
    }

    @Test
    public void testTrailingExtendCleanMatch()
    {
        final TestGenome genome = genome().set(CHR1, 131, 10, 'C');

        final SoftclipTailExtender ext = extender(genome);
        final TailExtensionResult res = ext.tryExtend(CHR1, 101, "30M10S", bases("A".repeat(30) + "C".repeat(10)));

        assertTrue(res.extended());
        assertEquals("40M", res.newCigar());
        assertEquals(101, res.newStart());
        assertEquals(10, res.basesExtendedTrail());
        assertEquals(0, res.basesExtendedLead());

        // statistics accumulate from the same call
        assertEquals(1, ext.statistics().recordsEvaluated());
        assertEquals(1, ext.statistics().recordsExtended());
        assertEquals(10, ext.statistics().basesExtendedTrail());
        assertEquals(0, ext.statistics().basesExtendedLead());
    }

    @Test
    public void testTrailingExtendWithOneMismatchInTen()
    {
        // one mismatch in 10 -> floor(10/10)=1 allowed -> accepts full 10
        final TestGenome genome = genome().set(CHR1, 131, 10, 'C').set(CHR1, 136, "G");   // single mismatch mid-tail

        final TailExtensionResult res = extender(genome)
                .tryExtend(CHR1, 101, "30M10S", bases("A".repeat(30) + "C".repeat(10)));

        assertTrue(res.extended());
        assertEquals("40M", res.newCigar());
        assertEquals(10, res.basesExtendedTrail());
    }

    @Test
    public void testTrailingExtendPartialThenMismatchBurst()
    {
        // 6 clean matches then a mismatch burst with no recovery; score peaks at 6 -> reclaim 6
        final TestGenome genome = genome().set(CHR1, 131, 6, 'C');

        final TailExtensionResult res = extender(genome)
                .tryExtend(CHR1, 101, "30M10S", bases("A".repeat(30) + "C".repeat(10)));

        assertTrue(res.extended());
        assertEquals("36M4S", res.newCigar());
        assertEquals(6, res.basesExtendedTrail());
    }

    @Test
    public void testTrailingExtendRejectShortAllMismatches()
    {
        // 3-base all-mismatch tail - below MinExtension=3, no extension
        final TailExtensionResult res = extender(genome())
                .tryExtend(CHR1, 101, "30M3S", bases("A".repeat(30) + "GGG"));

        assertFalse(res.extended());
    }

    @Test
    public void testLeadingExtendCleanMatch()
    {
        final TestGenome genome = genome().set(CHR1, 101, 10, 'T');

        final TailExtensionResult res = extender(genome)
                .tryExtend(CHR1, 111, "10S30M", bases("T".repeat(10) + "A".repeat(30)));

        assertTrue(res.extended());
        assertEquals("40M", res.newCigar());
        assertEquals(101, res.newStart());
        assertEquals(10, res.basesExtendedLead());
        assertEquals(0, res.basesExtendedTrail());
    }

    @Test
    public void testLeadingExtendPartial()
    {
        // last 5 of 10 leading softclip bases match, then mismatches with no recovery; reclaim 5
        final TestGenome genome = genome().set(CHR1, 106, 5, 'G');

        final TailExtensionResult res = extender(genome)
                .tryExtend(CHR1, 111, "10S30M", bases("C".repeat(5) + "G".repeat(5) + "A".repeat(30)));

        assertTrue(res.extended());
        assertEquals("5S35M", res.newCigar());
        assertEquals(106, res.newStart());
        assertEquals(5, res.basesExtendedLead());
    }

    @Test
    public void testBothEndsInOneCall()
    {
        final TestGenome genome = genome().set(CHR1, 106, 5, 'T').set(CHR1, 141, 5, 'C');

        final TailExtensionResult res = extender(genome)
                .tryExtend(CHR1, 111, "5S30M5S", bases("T".repeat(5) + "A".repeat(30) + "C".repeat(5)));

        assertTrue(res.extended());
        assertEquals("40M", res.newCigar());
        assertEquals(106, res.newStart());
        assertEquals(5, res.basesExtendedLead());
        assertEquals(5, res.basesExtendedTrail());
    }

    @Test
    public void testExtendsAcrossNCigar()
    {
        final TestGenome genome = new TestGenome().with(CHR1, 500, 'A').set(CHR1, 301, 5, 'C');

        final TailExtensionResult res = extender(genome)
                .tryExtend(CHR1, 101, "50M100N50M5S", bases("A".repeat(100) + "C".repeat(5)));

        assertTrue(res.extended());
        assertEquals("50M100N55M", res.newCigar());
        assertEquals(5, res.basesExtendedTrail());
    }

    @Test
    public void testMaxExtensionCap()
    {
        // MaxExtension=30 caps even when all 60 bases match
        final TestGenome genome = new TestGenome().with(CHR1, 500, 'A').set(CHR1, 131, 70, 'C');

        final TailExtensionResult res = extender(genome)
                .tryExtend(CHR1, 101, "30M60S", bases("A".repeat(30) + "C".repeat(60)));

        assertTrue(res.extended());
        assertEquals("60M30S", res.newCigar());
        assertEquals(30, res.basesExtendedTrail());
    }

    @Test
    public void testJunctionGuardDefersWhenGenomicDivergesTrailing()
    {
        // clip diverges from genomic ref at the annotated intron start -> defer to junction rescue
        final TestGenome genome = genome().set(CHR1, 131, 3, 'C');
        final Set<ChrBaseRegion> introns = new HashSet<>(Collections.singletonList(
                new ChrBaseRegion(CHR1, 131, 200)));

        final SoftclipTailExtender ext = extender(genome, introns);
        final TailExtensionResult res = ext.tryExtend(CHR1, 101, "30M10S", bases("A".repeat(30) + "C".repeat(10)));

        assertFalse(res.extended());
        assertEquals(1, ext.statistics().skippedForJunctionGuard());
    }

    @Test
    public void testRetainedIntronExtendsThroughJunctionTrailing()
    {
        // clip matches contiguous genome across the annotated intron - intron retention, extend not defer
        final TestGenome genome = genome().set(CHR1, 131, 10, 'C');
        final Set<ChrBaseRegion> introns = new HashSet<>(Collections.singletonList(
                new ChrBaseRegion(CHR1, 131, 200)));

        final SoftclipTailExtender ext = extender(genome, introns);
        final TailExtensionResult res = ext.tryExtend(CHR1, 101, "30M10S", bases("A".repeat(30) + "C".repeat(10)));

        assertTrue(res.extended());
        assertEquals("40M", res.newCigar());
        assertEquals(0, ext.statistics().skippedForJunctionGuard());
    }

    @Test
    public void testJunctionGuardDefersWhenGenomicDivergesLeading()
    {
        // clip diverges from genomic ref at the annotated intron end -> defer
        final TestGenome genome = genome().set(CHR1, 108, 3, 'T');
        final Set<ChrBaseRegion> introns = new HashSet<>(Collections.singletonList(
                new ChrBaseRegion(CHR1, 21, 110)));

        final SoftclipTailExtender ext = extender(genome, introns);
        final TailExtensionResult res = ext.tryExtend(CHR1, 111, "10S30M", bases("T".repeat(10) + "A".repeat(30)));

        assertFalse(res.extended());
        assertEquals(1, ext.statistics().skippedForJunctionGuard());
    }

    @Test
    public void testRetainedIntronExtendsThroughJunctionLeading()
    {
        // clip cleanly matches genome across the annotated intron end - retention, extend not defer
        final TestGenome genome = genome().set(CHR1, 101, 10, 'T');
        final Set<ChrBaseRegion> introns = new HashSet<>(Collections.singletonList(
                new ChrBaseRegion(CHR1, 21, 110)));

        final SoftclipTailExtender ext = extender(genome, introns);
        final TailExtensionResult res = ext.tryExtend(CHR1, 111, "10S30M", bases("T".repeat(10) + "A".repeat(30)));

        assertTrue(res.extended());
        assertEquals("40M", res.newCigar());
        assertEquals(101, res.newStart());
        assertEquals(0, ext.statistics().skippedForJunctionGuard());
    }

    @Test
    public void testNoOpAndRefusalGuards()
    {
        // no softclip -> unchanged
        assertFalse(extender(genome()).tryExtend(CHR1, 101, "30M", repeatedBase(30, 'A')).extended());

        // 2S < MinSoftclipLength=3 -> unchanged
        assertFalse(extender(genome().set(CHR1, 131, 2, 'C'))
                .tryExtend(CHR1, 101, "30M2S", bases("A".repeat(30) + "CC")).extended());

        // hard clip -> skipped as complex shape
        final SoftclipTailExtender hardClip = extender(genome());
        assertFalse(hardClip.tryExtend(CHR1, 101, "10H30M5S", repeatedBase(35, 'A')).extended());
        assertEquals(1, hardClip.statistics().skippedComplexShape());

        // op adjacent to the trailing S is I, not M -> refuse as complex shape
        final SoftclipTailExtender indelAdjacent = extender(genome().set(CHR1, 131, 10, 'C'));
        assertFalse(indelAdjacent.tryExtend(CHR1, 101, "30M5I3S", bases("A".repeat(35) + "CCC")).extended());
        assertTrue(indelAdjacent.statistics().skippedComplexShape() >= 1);

        // disabled config -> no-op even with a clean match
        final SoftclipTailExtender disabled = new SoftclipTailExtender(
                genome().set(CHR1, 131, 10, 'C').asRefSource(), null, TailExtensionConfig.defaults());
        assertFalse(disabled.tryExtend(CHR1, 101, "30M10S", bases("A".repeat(30) + "C".repeat(10))).extended());

        // null chromosome / cigar / read bases -> no-op
        final SoftclipTailExtender nullInputs = extender(new TestGenome());
        assertFalse(nullInputs.tryExtend(null, 101, "30M5S", repeatedBase(35, 'A')).extended());
        assertFalse(nullInputs.tryExtend(CHR1, 101, null, repeatedBase(35, 'A')).extended());
        assertFalse(nullInputs.tryExtend(CHR1, 101, "30M5S", null).extended());
    }
}
