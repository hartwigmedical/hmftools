package com.hartwig.hmftools.tars.liftback.tailextend;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import com.hartwig.hmftools.tars.common.SpliceCommon;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;

import org.junit.Test;

public class TerminalMicroJunctionCollapserTest
{
    private static final String CHR1 = "chr1";

    // single-chromosome in-memory ref backed by the shared fixture (MockRefGenome, 1-based inclusive).
    private static RefSequenceSource refSource(final byte[] chr1)
    {
        return TarsTestFixtures.refSource(CHR1, new String(chr1, StandardCharsets.US_ASCII));
    }

    private static TerminalMicroJunctionCollapser collapser(final byte[] chr1)
    {
        return new TerminalMicroJunctionCollapser(refSource(chr1), SpliceCommon.MIN_JUNCTION_ANCHOR);
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
        // 1M10N8M: leading 1bp anchor continues the near exon backwards (pos 111); near exon is a full
        // 8bp block (>= threshold) so only the leading anchor collapses -> 9M, start moved back to 111.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        chr1[100] = 'C';                                       // pos 101: spurious far base
        // pos 111 = contiguous anchor base, 112..119 = the 8bp near exon
        final byte[] window = { 'C', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T' };  // 111..119
        System.arraycopy(window, 0, chr1, 110, window.length);
        final byte[] read = { 'C', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T' };    // anchor + 8M near exon

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "1M10N8M", read);

        assertTrue(res.Collapsed);
        assertEquals("9M", res.NewCigar);
        assertEquals(111, res.NewStart);
    }

    @Test
    public void testSubAnchorJunctionClipped()
    {
        // 2M80N149M where the 2bp anchor matches neither the far exon nor the contiguous near exon (split and
        // contiguous both ~0). Contiguous ties, so the sub-threshold intron is dropped and the anchor clipped
        // -> 2S149M, start moved to the near exon.
        final byte[] chr1 = new byte[400];
        Arrays.fill(chr1, (byte) 'T');
        chr1[80] = 'G'; chr1[81] = 'G';     // pos 81,82: intron bases, do not match read's "AC"
        final byte[] read = new byte[151];
        Arrays.fill(read, (byte) 'A');
        read[0] = 'A'; read[1] = 'C';

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 1, "2M80N149M", read);

        assertTrue(res.Collapsed);
        assertEquals("2S149M", res.NewCigar);
        assertEquals(83, res.NewStart);   // nearStart = 1 + 2 + 80
    }

    @Test
    public void testRealShortJunctionKeptWhenSplitScoresBetter()
    {
        // 50M200N5M: the 5bp anchor matches the far exon (bwa's placement, 351..355) exactly but does NOT
        // continue the near exon contiguously (151..155 mismatch). Split (+5) beats contiguous (0), so it is a
        // real short-overhang junction and is kept.
        final byte[] chr1 = new byte[400];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] anchor = { 'A', 'C', 'G', 'A', 'C' };
        System.arraycopy(anchor, 0, chr1, 350, anchor.length);   // far exon 351..355 matches the anchor
        // near-contiguous 151..155 stay 'T' -> anchor doesn't extend contiguously
        final byte[] read = new byte[55];
        Arrays.fill(read, (byte) 'G');
        System.arraycopy(anchor, 0, read, 50, anchor.length);

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "50M200N5M", read);

        assertFalse(res.Collapsed);
    }

    @Test
    public void testCleanFarAnchorWithJunkClipKept()
    {
        // exp8 shape: 50M200N5M3S. The 5M matches the far exon cleanly (split +5); the 3S matches neither the
        // far continuation nor the contiguous intron (contiguous 0). Split wins -> junction kept, clip stays.
        final byte[] chr1 = new byte[400];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] anchor = { 'A', 'C', 'G', 'A', 'C' };
        System.arraycopy(anchor, 0, chr1, 350, anchor.length);   // far exon 351..355 matches the 5M anchor
        final byte[] read = new byte[58];
        Arrays.fill(read, (byte) 'G');
        System.arraycopy(anchor, 0, read, 50, anchor.length);    // 5M anchor; last 3 ('G') = junk clip

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "50M200N5M3S", read);

        assertFalse(res.Collapsed);
    }

    @Test
    public void testAnchorAtThresholdKept()
    {
        // both flanks meet MIN_JUNCTION_ANCHOR (8): a trusted junction, neither end is a collapse candidate.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'C');
        final byte[] read = new byte[16];
        Arrays.fill(read, (byte) 'C');

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "8M10N8M", read);

        assertFalse(res.Collapsed);
    }

    @Test
    public void testTrailingNonContiguousClipped()
    {
        // Anchor base doesn't continue the near exon (pos 106 is 'A', read has 'C'): reclaim nothing, drop
        // the intron, clip the anchor -> 5M1S.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'A');
        chr1[116] = 'C';   // far base at 117 matches; contiguous pos 106 stays 'A'
        final byte[] read = { 'A', 'A', 'A', 'A', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M", read);

        assertTrue(res.Collapsed);
        assertEquals("5M1S", res.NewCigar);
        assertEquals(101, res.NewStart);
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
        // Anchor + first 3 softclip bases continue the near exon, last 2 diverge with no recovery; score
        // peaks at 3 reclaimed past the 5M near exon -> 8M2S.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] genome = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'G', 'G' };   // 109,110 diverge
        System.arraycopy(genome, 0, chr1, 100, genome.length);
        final byte[] read = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M4S", read);

        assertTrue(res.Collapsed);
        assertEquals("8M2S", res.NewCigar);
        assertEquals(101, res.NewStart);
    }

    @Test
    public void testTrailingSoftclipAnchorMismatchClipped()
    {
        // Anchor base mismatches the near-exon continuation (pos 106 = 'T'): score never goes positive, so
        // the whole terminal region is clipped and the intron dropped -> 5M5S.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] nearExon = { 'A', 'C', 'G', 'T', 'A' };
        System.arraycopy(nearExon, 0, chr1, 100, nearExon.length);   // 101..105, 106 stays 'T'
        final byte[] read = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "5M10N1M4S", read);

        assertTrue(res.Collapsed);
        assertEquals("5M5S", res.NewCigar);
        assertEquals(101, res.NewStart);
    }

    @Test
    public void testLeadingSoftclipReclaimed()
    {
        // 4S1M10N8M: anchor + 4 softclip bases continue the 8bp near exon backwards contiguously; the near
        // exon (>= threshold) isn't itself a collapse candidate, so only the leading tail folds -> 13M, start 107.
        final byte[] chr1 = new byte[200];
        Arrays.fill(chr1, (byte) 'T');
        final byte[] window = { 'C', 'G', 'T', 'A', 'C' };       // 107..111: softclip(4) + anchor(1)
        System.arraycopy(window, 0, chr1, 106, window.length);
        final byte[] nearExon = { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T' };   // 112..119, 8bp
        System.arraycopy(nearExon, 0, chr1, 111, nearExon.length);
        final byte[] read = { 'C', 'G', 'T', 'A', 'C', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T' };

        final TerminalCollapseResult res = collapser(chr1).tryCollapse(CHR1, 101, "4S1M10N8M", read);

        assertTrue(res.Collapsed);
        assertEquals("13M", res.NewCigar);
        assertEquals(107, res.NewStart);
    }
}
