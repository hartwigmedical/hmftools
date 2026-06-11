package com.hartwig.hmftools.redux.splice.tailextend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import com.hartwig.hmftools.redux.splice.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.redux.splice.rescue.ChrIntron;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;

public class SoftclipTailExtenderTest
{
    private static final String CHR1 = "chr1";

    private static RefSequenceSource refSource(final Map<String, byte[]> ref)
    {
        return (chrom, posStart, posEnd) ->
        {
            final byte[] bases = ref.get(chrom);
            if(bases == null)
                return null;
            final int startIdx = posStart - 1;
            final int endIdx = posEnd - 1;
            if(startIdx < 0 || endIdx >= bases.length || endIdx < startIdx)
                return null;
            return Arrays.copyOfRange(bases, startIdx, endIdx + 1);
        };
    }

    private static byte[] fill(final int length, final char base)
    {
        final byte[] out = new byte[length];
        Arrays.fill(out, (byte) base);
        return out;
    }

    private static SoftclipTailExtender extender(final Map<String, byte[]> ref, final Set<ChrIntron> introns)
    {
        return new SoftclipTailExtender(
                refSource(ref),
                introns == null ? null : new AnnotatedJunctionIndex(introns),
                TailExtensionConfig.enabledDefaults());
    }

    private static SoftclipTailExtender extender(final Map<String, byte[]> ref)
    {
        return extender(ref, null);
    }

    @Test
    public void testTrailingExtendCleanMatch()
    {
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 140; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[40];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "30M10S", readBases);

        assertTrue(res.Extended);
        assertEquals("40M", res.NewCigar);
        assertEquals(101, res.NewStart);
        assertEquals(10, res.BasesExtendedTrail);
        assertEquals(0, res.BasesExtendedLead);
    }

    @Test
    public void testTrailingExtendWithOneMismatchInTen()
    {
        // one mismatch in 10 -> floor(10/10)=1 allowed -> accepts full 10
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 140; ++i) chr1[i] = 'C';
        chr1[135] = 'G';                                 // single mismatch in the middle of the tail
        final byte[] readBases = new byte[40];
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "30M10S", readBases);

        assertTrue(res.Extended);
        assertEquals("40M", res.NewCigar);
        assertEquals(10, res.BasesExtendedTrail);
    }

    @Test
    public void testTrailingExtendPartialThenMismatchBurst()
    {
        // 6 clean matches then a mismatch burst with no recovery; score peaks at 6 -> reclaim 6
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 136; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[40];
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "30M10S", readBases);

        assertTrue(res.Extended);
        assertEquals("36M4S", res.NewCigar);
        assertEquals(6, res.BasesExtendedTrail);
    }

    @Test
    public void testTrailingExtendRejectShortAllMismatches()
    {
        // 3-base all-mismatch tail — below MinExtension=3, no extension
        final byte[] chr1 = fill(200, 'A');
        final byte[] readBases = new byte[33];
        for(int i = 0; i < 30; ++i) readBases[i] = 'A';
        for(int i = 30; i < 33; ++i) readBases[i] = 'G';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "30M3S", readBases);

        assertFalse(res.Extended);
    }

    @Test
    public void testLeadingExtendCleanMatch()
    {
        final byte[] chr1 = fill(200, 'A');
        for(int i = 100; i < 110; ++i) chr1[i] = 'T';
        final byte[] readBases = new byte[40];
        for(int i = 0; i < 10; ++i) readBases[i] = 'T';
        for(int i = 10; i < 40; ++i) readBases[i] = 'A';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 111, "10S30M", readBases);

        assertTrue(res.Extended);
        assertEquals("40M", res.NewCigar);
        assertEquals(101, res.NewStart);
        assertEquals(10, res.BasesExtendedLead);
        assertEquals(0, res.BasesExtendedTrail);
    }

    @Test
    public void testLeadingExtendPartial()
    {
        // last 5 of 10 leading softclip bases match, then mismatches with no recovery; reclaim 5
        final byte[] chr1 = fill(200, 'A');
        for(int i = 105; i < 110; ++i) chr1[i] = 'G';
        final byte[] readBases = new byte[40];
        for(int i = 0; i < 5; ++i) readBases[i] = 'C';
        for(int i = 5; i < 10; ++i) readBases[i] = 'G';
        for(int i = 10; i < 40; ++i) readBases[i] = 'A';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 111, "10S30M", readBases);

        assertTrue(res.Extended);
        assertEquals("5S35M", res.NewCigar);
        assertEquals(106, res.NewStart);
        assertEquals(5, res.BasesExtendedLead);
    }

    @Test
    public void testBothEndsInOneCall()
    {
        final byte[] chr1 = fill(200, 'A');
        for(int i = 105; i < 110; ++i) chr1[i] = 'T';
        for(int i = 140; i < 145; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[40];
        for(int i = 0; i < 5; ++i) readBases[i] = 'T';
        for(int i = 5; i < 35; ++i) readBases[i] = 'A';
        for(int i = 35; i < 40; ++i) readBases[i] = 'C';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 111, "5S30M5S", readBases);

        assertTrue(res.Extended);
        assertEquals("40M", res.NewCigar);
        assertEquals(106, res.NewStart);
        assertEquals(5, res.BasesExtendedLead);
        assertEquals(5, res.BasesExtendedTrail);
    }

    @Test
    public void testExtendsAcrossNCigar()
    {
        final byte[] chr1 = fill(500, 'A');
        for(int i = 300; i < 305; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[105];
        for(int i = 100; i < 105; ++i) readBases[i] = 'C';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "50M100N50M5S", readBases);

        assertTrue(res.Extended);
        assertEquals("50M100N55M", res.NewCigar);
        assertEquals(5, res.BasesExtendedTrail);
    }

    @Test
    public void testMaxExtensionCap()
    {
        // MaxExtension=30 caps even when all 60 bases match
        final byte[] chr1 = fill(500, 'A');
        for(int i = 130; i < 200; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[90];
        for(int i = 30; i < 90; ++i) readBases[i] = 'C';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "30M60S", readBases);

        assertTrue(res.Extended);
        assertEquals("60M30S", res.NewCigar);
        assertEquals(30, res.BasesExtendedTrail);
    }

    @Test
    public void testJunctionGuardDefersWhenGenomicDivergesTrailing()
    {
        // clip diverges from genomic ref at the annotated intron start -> defer to junction rescue
        final byte[] chr1 = fill(200, 'A');
        chr1[130] = 'C'; chr1[131] = 'C'; chr1[132] = 'C';
        final byte[] readBases = new byte[40];
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';

        final Set<ChrIntron> introns = new HashSet<>(Collections.singletonList(
                new ChrIntron(CHR1, 131, 200)));

        final SoftclipTailExtender ext = extender(Collections.singletonMap(CHR1, chr1), introns);
        final TailExtensionResult res = ext.tryExtend(CHR1, 101, "30M10S", readBases);

        assertFalse(res.Extended);
        assertEquals(1, ext.statistics().skippedForJunctionGuard());
    }

    @Test
    public void testRetainedIntronExtendsThroughJunctionTrailing()
    {
        // clip matches contiguous genome across the annotated intron — intron retention, extend not defer
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 140; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[40];
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';

        final Set<ChrIntron> introns = new HashSet<>(Collections.singletonList(
                new ChrIntron(CHR1, 131, 200)));

        final SoftclipTailExtender ext = extender(Collections.singletonMap(CHR1, chr1), introns);
        final TailExtensionResult res = ext.tryExtend(CHR1, 101, "30M10S", readBases);

        assertTrue(res.Extended);
        assertEquals("40M", res.NewCigar);
        assertEquals(0, ext.statistics().skippedForJunctionGuard());
    }

    @Test
    public void testJunctionGuardDefersWhenGenomicDivergesLeading()
    {
        // clip diverges from genomic ref at the annotated intron end -> defer
        final byte[] chr1 = fill(200, 'A');
        chr1[107] = 'T'; chr1[108] = 'T'; chr1[109] = 'T';
        final byte[] readBases = new byte[40];
        for(int i = 0; i < 10; ++i) readBases[i] = 'T';

        final Set<ChrIntron> introns = new HashSet<>(Collections.singletonList(
                new ChrIntron(CHR1, 21, 110)));

        final SoftclipTailExtender ext = extender(Collections.singletonMap(CHR1, chr1), introns);
        final TailExtensionResult res = ext.tryExtend(CHR1, 111, "10S30M", readBases);

        assertFalse(res.Extended);
        assertEquals(1, ext.statistics().skippedForJunctionGuard());
    }

    @Test
    public void testRetainedIntronExtendsThroughJunctionLeading()
    {
        // clip cleanly matches genome across the annotated intron end — retention, extend not defer
        final byte[] chr1 = fill(200, 'A');
        for(int i = 100; i < 110; ++i) chr1[i] = 'T';
        final byte[] readBases = new byte[40];
        for(int i = 0; i < 10; ++i) readBases[i] = 'T';

        final Set<ChrIntron> introns = new HashSet<>(Collections.singletonList(
                new ChrIntron(CHR1, 21, 110)));

        final SoftclipTailExtender ext = extender(Collections.singletonMap(CHR1, chr1), introns);
        final TailExtensionResult res = ext.tryExtend(CHR1, 111, "10S30M", readBases);

        assertTrue(res.Extended);
        assertEquals("40M", res.NewCigar);
        assertEquals(101, res.NewStart);
        assertEquals(0, ext.statistics().skippedForJunctionGuard());
    }

    @Test
    public void testNoSoftclipUnchanged()
    {
        final byte[] chr1 = fill(200, 'A');
        final byte[] readBases = fill(30, 'A');

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "30M", readBases);

        assertFalse(res.Extended);
    }

    @Test
    public void testShortSoftclipBelowMinUnchanged()
    {
        // 2S < MinSoftclipLength=3 -> unchanged
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 132; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[32];
        readBases[30] = 'C'; readBases[31] = 'C';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "30M2S", readBases);

        assertFalse(res.Extended);
    }

    @Test
    public void testHardClipRefusal()
    {
        // hard clip -> skipped as complex shape
        final byte[] chr1 = fill(200, 'A');
        final byte[] readBases = fill(35, 'A');

        final SoftclipTailExtender ext = extender(Collections.singletonMap(CHR1, chr1));
        final TailExtensionResult res = ext.tryExtend(CHR1, 101, "10H30M5S", readBases);

        assertFalse(res.Extended);
        assertEquals(1, ext.statistics().skippedComplexShape());
    }

    @Test
    public void testIndelAdjacentToSoftclipRefused()
    {
        // op adjacent to trailing S is I, not M -> refuse
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 140; ++i) chr1[i] = 'C';
        final byte[] readBases = fill(38, 'A');
        for(int i = 35; i < 38; ++i) readBases[i] = 'C';

        final SoftclipTailExtender ext = extender(Collections.singletonMap(CHR1, chr1));
        final TailExtensionResult res = ext.tryExtend(CHR1, 101, "30M5I3S", readBases);

        assertFalse(res.Extended);
        assertTrue(ext.statistics().skippedComplexShape() >= 1);
    }

    @Test
    public void testDisabledConfigIsNoOp()
    {
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 140; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[40];
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';

        final SoftclipTailExtender ext = new SoftclipTailExtender(
                refSource(Collections.singletonMap(CHR1, chr1)), null, TailExtensionConfig.defaults());
        final TailExtensionResult res = ext.tryExtend(CHR1, 101, "30M10S", readBases);

        assertFalse(res.Extended);
    }

    @Test
    public void testStatisticsCounters()
    {
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 140; ++i) chr1[i] = 'C';
        final byte[] readBases = new byte[40];
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';

        final SoftclipTailExtender ext = extender(Collections.singletonMap(CHR1, chr1));
        ext.tryExtend(CHR1, 101, "30M10S", readBases);
        ext.tryExtend(CHR1, 101, "40M", readBases);     // no softclip — still counted as evaluated

        assertEquals(2, ext.statistics().recordsEvaluated());
        assertEquals(1, ext.statistics().recordsExtended());
        assertEquals(10, ext.statistics().basesExtendedTrail());
        assertEquals(0, ext.statistics().basesExtendedLead());
    }

    @Test
    public void testNullInputsAreNoOp()
    {
        final SoftclipTailExtender ext = extender(new HashMap<>());
        assertFalse(ext.tryExtend(null, 101, "30M5S", fill(35, 'A')).Extended);
        assertFalse(ext.tryExtend(CHR1, 101, null, fill(35, 'A')).Extended);
        assertFalse(ext.tryExtend(CHR1, 101, "30M5S", null).Extended);
    }
}
