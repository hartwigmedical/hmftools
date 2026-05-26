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

    // Builds a RefSequenceSource backed by an in-memory map. Coords are 1-based inclusive. Returns
    // null when the chromosome is unknown or the range falls outside the array.
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
        // Primary 30M10S at chr1:101. Trailing 10 read bases are "CCCCCCCCCC". Ref at chr1:131..140
        // = "CCCCCCCCCC". Should fully extend the trailing 10S into M.
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
        // 10 trailing bases, one mismatch -> floor(10/10)=1 allowed -> still accepts full 10.
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
        // First 6 bases match, then a burst of mismatches. The max(1, length/10) floor allows the
        // 7th base to be a mismatch (m=1, allowed=1) -> bestLength=7. The 8th and 9th positions
        // would push m to 2 and 3 which exceed allowed+1 -> break.
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 136; ++i) chr1[i] = 'C';     // ref C C C C C C ...
        // remaining positions stay A; read has C there -> mismatches at positions 6+ of softclip
        final byte[] readBases = new byte[40];
        for(int i = 30; i < 40; ++i) readBases[i] = 'C';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 101, "30M10S", readBases);

        assertTrue(res.Extended);
        assertEquals("37M3S", res.NewCigar);
        assertEquals(7, res.BasesExtendedTrail);
    }

    @Test
    public void testTrailingExtendRejectShortAllMismatches()
    {
        // 3-base tail, every base mismatches. Best matched run is 1 (covered by the floor-of-1
        // tolerance), which is below MinExtension=3 → no extension.
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
        // Primary 10S30M at chr1:111. Leading 10 read bases = "TTTTTTTTTT". Ref at chr1:101..110
        // = "TTTTTTTTTT". Should extend, new start = 101.
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
        // Leading 10S, only the last 5 bases of the softclip match ref. The walk-back algorithm
        // pairs read[softclipLen-1-i] with ref[refLen-1-i]. After 5 clean matches, the 6th base
        // is the first mismatch (m=1, allowed=1, accept) -> bestLength=6. The 7th pushes m=2
        // which exceeds allowed+1, so it breaks.
        final byte[] chr1 = fill(200, 'A');
        for(int i = 105; i < 110; ++i) chr1[i] = 'G';     // ref[106..110] = G
        final byte[] readBases = new byte[40];
        for(int i = 0; i < 5; ++i) readBases[i] = 'C';     // first 5 softclip bases mismatch
        for(int i = 5; i < 10; ++i) readBases[i] = 'G';    // last 5 match
        for(int i = 10; i < 40; ++i) readBases[i] = 'A';

        final TailExtensionResult res = extender(Collections.singletonMap(CHR1, chr1))
                .tryExtend(CHR1, 111, "10S30M", readBases);

        assertTrue(res.Extended);
        assertEquals("4S36M", res.NewCigar);
        assertEquals(105, res.NewStart);
        assertEquals(6, res.BasesExtendedLead);
    }

    @Test
    public void testBothEndsInOneCall()
    {
        // Lead 5S + trail 5S, both with ref support.
        final byte[] chr1 = fill(200, 'A');
        for(int i = 105; i < 110; ++i) chr1[i] = 'T';     // ref for leading (pos 106..110)
        for(int i = 140; i < 145; ++i) chr1[i] = 'C';     // ref for trailing (pos 141..145)
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
        // 50M100N50M5S — internal N op, trailing 5S. Extension touches only the trailing side.
        // Ref span = 50 + 100 + 50 = 200 -> alignmentEnd = 101 + 200 - 1 = 300.
        // Ref at 301..305 = "CCCCC", softclip = "CCCCC".
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
        // 60S trailing softclip, ref matches all 60 bases. MaxExtension=30 caps extension at 30.
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
    public void testJunctionGuardBlocksTrailing()
    {
        // Annotated intron starts at chr1:131 (just past the M run end). Even though ref matches,
        // the guard must skip — this is L1 rescue territory.
        final byte[] chr1 = fill(200, 'A');
        for(int i = 130; i < 140; ++i) chr1[i] = 'C';
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
    public void testJunctionGuardBlocksLeading()
    {
        // Annotated intron ends at chr1:110 (just before the M run start). Guard must skip.
        final byte[] chr1 = fill(200, 'A');
        for(int i = 100; i < 110; ++i) chr1[i] = 'T';
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
        // 2S softclip < MinSoftclipLength=3 -> unchanged.
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
        // Hard clip present -> skip with countSkippedComplexShape.
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
        // 30M5I3S — the op adjacent to the trailing S is I, not M. Refuse to extend.
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
        ext.tryExtend(CHR1, 101, "40M", readBases);     // no softclip, eval still counted

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
