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
        // Read aligned 5M, then a 10N intron, then a 1bp terminal anchor (5M10N1M). The 6th read base
        // also matches the genome continuing the first exon (pos 106), so the junction is unnecessary
        // and collapses to 6M. Ref: exon1 at 101..106 = "ACGTAC" (the 6th base C continues it), intron
        // 107..116, the spurious far base at 117 also 'C'.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        // 101..106 = ACGTAC  (read M run + the contiguous 6th base)
        final byte[] exon = { 'A', 'C', 'G', 'T', 'A', 'C' };
        System.arraycopy(exon, 0, chr1, 100, 6);
        chr1[116] = 'C';   // pos 117: the base the contig walk over-extended onto across the intron

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
        // Mirror on the leading end: 1M 10N 5M. The 1 leading base also matches the genome immediately
        // before the second exon (contiguous), so collapse to 6M with the start moved back by 1+10.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        // near exon (the 5M) sits at 112..116; the contiguous base before it (111) matches the read's
        // first base. far spurious base at 101.
        chr1[100] = 'C';                                   // pos 101: spurious leading 1M
        final byte[] nearExon = { 'C', 'A', 'C', 'G', 'T' }; // pos 111..115: contiguous base + 5M start
        System.arraycopy(nearExon, 0, chr1, 110, 5);
        // read first base C matches both pos 101 (far) and pos 111 (contiguous continuation)
        final byte[] read = { 'C', 'A', 'C', 'G', 'T' };

        // 1M at 101, 10N (102..111), 4M at 112. nearStart = 101 + 1 + 10 = 112; check pos 111.
        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "1M10N4M", read);

        assertTrue(res.Collapsed);
        assertEquals("5M", res.NewCigar);
        assertEquals(111, res.NewStart);   // nearStart(112) - y(1)
    }

    @Test
    public void testLegitShortAnchorJunctionKept()
    {
        // The legit "2M 80N 149M" case: the read genuinely starts 2 bases into exon1. Those 2 leading
        // bases do NOT match the genome immediately before exon2 (an intron sits there), so the ref
        // check fails and the junction is KEPT. Anchor length alone never drives the decision.
        final byte[] chr1 = new byte[400];
        Arrays.fill(chr1, (byte) 'T');
        // read's first 2 bases = "AC" (real exon1 tail). The 2 bases before the near exon (the intron's
        // last 2 bases, pos 81..82 given start 1, 2M then 80N then exon at 83) are "GG" - no match.
        chr1[80] = 'G'; chr1[81] = 'G';     // pos 81,82: intron bases before exon2, do not match "AC"
        final byte[] read = new byte[151];
        Arrays.fill(read, (byte) 'A');
        read[0] = 'A'; read[1] = 'C';

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 1, "2M80N149M", read);

        assertFalse(res.Collapsed);
    }

    @Test
    public void testAnchorAboveMaxNotCollapsed()
    {
        // A 3bp terminal anchor is above the max (2), so it is never a collapse candidate regardless of
        // whether the bases happen to match contiguously.
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
        // Terminal 1M anchor, but the over-extended base does NOT match the contiguous continuation of
        // the near exon (pos 106 is 'A', read's 6th base is 'C') - so the junction has real support and
        // is kept.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'A');
        chr1[116] = 'C';   // far base at 117 matches; contiguous pos 106 stays 'A'
        final byte[] read = { 'A', 'A', 'A', 'A', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M", read);

        assertFalse(res.Collapsed);
    }
}
