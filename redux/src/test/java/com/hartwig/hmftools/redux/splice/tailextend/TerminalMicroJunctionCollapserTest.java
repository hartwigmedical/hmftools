package com.hartwig.hmftools.redux.splice.tailextend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Map;

import org.junit.Test;

import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;

public class TerminalMicroJunctionCollapserTest
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

    private static TerminalMicroJunctionCollapser collapser(final byte[] chr1)
    {
        return new TerminalMicroJunctionCollapser(refSource(Map.of(CHR1, chr1)), 2);
    }

    @Test
    public void testTrailingSpuriousMicroJunctionCollapsed()
    {
        // 5M10N1M: anchor base matches the contiguous continuation of the near exon (pos 106), so collapses to 6M.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        // 101..106 = ACGTAC
        final byte[] exon = { 'A', 'C', 'G', 'T', 'A', 'C' };
        System.arraycopy(exon, 0, chr1, 100, 6);
        chr1[116] = 'C';   // pos 117: spurious far base

        // read bases: first 5 = ACGTA (5M), 6th = C (the over-extended 1M, also matches pos 106)
        final byte[] read = { 'A', 'C', 'G', 'T', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M", read);

        assertTrue(res.Collapsed);
        assertEquals("6M", res.NewCigar);
        assertEquals(101, res.NewStart);
    }

    @Test
    public void testLeadingSpuriousMicroJunctionCollapsed()
    {
        // 1M10N4M: anchor base matches contiguous genome before the near exon; collapses to 5M, start -11.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        chr1[100] = 'C';                                       // pos 101: spurious far base
        final byte[] nearExon = { 'C', 'A', 'C', 'G', 'T' }; // pos 111..115
        System.arraycopy(nearExon, 0, chr1, 110, 5);
        final byte[] read = { 'C', 'A', 'C', 'G', 'T' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "1M10N4M", read);

        assertTrue(res.Collapsed);
        assertEquals("5M", res.NewCigar);
        assertEquals(111, res.NewStart);
    }

    @Test
    public void testLegitShortAnchorJunctionKept()
    {
        // Genuine 2M80N149M: anchor bases don't match the contiguous genome before exon2 (intron sits there).
        final byte[] chr1 = new byte[400];
        Arrays.fill(chr1, (byte) 'T');
        chr1[80] = 'G'; chr1[81] = 'G';     // pos 81,82: intron bases, do not match read's "AC"
        final byte[] read = new byte[151];
        Arrays.fill(read, (byte) 'A');
        read[0] = 'A'; read[1] = 'C';

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 1, "2M80N149M", read);

        assertFalse(res.Collapsed);
    }

    @Test
    public void testAnchorAboveMaxNotCollapsed()
    {
        // 3bp anchor exceeds max (2), so never a collapse candidate.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'C');
        final byte[] read = new byte[9];
        Arrays.fill(read, (byte) 'C');

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N3M", read);

        assertFalse(res.Collapsed);
    }

    @Test
    public void testTrailingNonContiguousKept()
    {
        // Anchor base doesn't match contiguous genome (pos 106 is 'A', read has 'C'): junction kept.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'A');
        chr1[116] = 'C';   // far base at 117 matches; contiguous pos 106 stays 'A'
        final byte[] read = { 'A', 'A', 'A', 'A', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M", read);

        assertFalse(res.Collapsed);
    }

    @Test
    public void testTrailingSoftclipFullyReclaimed()
    {
        // 5M10N1M4S: anchor + all 4 softclip bases match contiguous genome → 10M.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] genome = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C' };
        System.arraycopy(genome, 0, chr1, 100, genome.length);   // 101..110
        final byte[] read = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M4S", read);

        assertTrue(res.Collapsed);
        assertEquals("10M", res.NewCigar);
        assertEquals(101, res.NewStart);
    }

    @Test
    public void testTrailingSoftclipPartiallyReclaimed()
    {
        // Last 2 softclip bases diverge beyond budget → anchor + 3 bases reclaimed (9M1S).
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] genome = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'G', 'G' };   // 109,110 diverge
        System.arraycopy(genome, 0, chr1, 100, genome.length);
        final byte[] read = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M4S", read);

        assertTrue(res.Collapsed);
        assertEquals("9M1S", res.NewCigar);
        assertEquals(101, res.NewStart);
    }

    @Test
    public void testTrailingSoftclipAnchorMismatchKept()
    {
        // Anchor base mismatches contiguous genome (pos 106 = 'T'): junction kept regardless of softclip.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] nearExon = { 'A', 'C', 'G', 'T', 'A' };
        System.arraycopy(nearExon, 0, chr1, 100, nearExon.length);   // 101..105, 106 stays 'T'
        final byte[] read = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M4S", read);

        assertFalse(res.Collapsed);
    }

    @Test
    public void testLeadingSoftclipReclaimed()
    {
        // 4S1M10N5M: anchor + 4 softclip bases match contiguous genome before near exon → 10M, start 107.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] window = { 'C', 'G', 'T', 'A', 'C' };       // 107..111: softclip(4) + anchor(1)
        System.arraycopy(window, 0, chr1, 106, window.length);
        final byte[] nearExon = { 'A', 'C', 'G', 'T', 'A' };     // 112..116
        System.arraycopy(nearExon, 0, chr1, 111, nearExon.length);
        final byte[] read = { 'C', 'G', 'T', 'A', 'C', 'A', 'C', 'G', 'T', 'A' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "4S1M10N5M", read);

        assertTrue(res.Collapsed);
        assertEquals("10M", res.NewCigar);
        assertEquals(107, res.NewStart);
    }
}
