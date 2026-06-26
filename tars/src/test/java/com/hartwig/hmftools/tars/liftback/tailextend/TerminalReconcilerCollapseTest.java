package com.hartwig.hmftools.tars.liftback.tailextend;

import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.bases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.tars.common.TarsConstants;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TestGenome;

import org.junit.Test;

public class TerminalReconcilerCollapseTest
{
    private static final String CHR1 = "chr1";

    private static TerminalReconciler collapser(final TestGenome genome)
    {
        return new TerminalReconciler(
                genome.asRefSource(), TarsConstants.MIN_JUNCTION_ANCHOR, null, TailExtensionConfig.defaults());
    }

    @Test
    public void testTrailingSpuriousMicroJunctionCollapsed()
    {
        // 5M10N1M: anchor base matches the contiguous continuation of the near exon (pos 106), so collapses to 6M.
        // 101..106 = ACGTAC; pos 117 = spurious far base.
        TestGenome genome = new TestGenome().with(CHR1, 200, 'T').set(CHR1, 101, "ACGTAC").set(CHR1, 117, "C");

        // read bases: first 5 = ACGTA (5M), 6th = C (the over-extended 1M, also matches pos 106)
        TerminalCollapseResult res = collapser(genome).tryCollapse(CHR1, 101, "5M10N1M", bases("ACGTAC"));

        assertTrue(res.collapsed());
        assertEquals("6M", res.newCigar());
        assertEquals(101, res.newStart());
    }

    @Test
    public void testLeadingSpuriousMicroJunctionCollapsed()
    {
        // 1M10N8M: leading 1bp anchor continues the near exon backwards (pos 111); near exon is a full
        // 8bp block (>= threshold) so only the leading anchor collapses -> 9M, start moved back to 111.
        // pos 101 = spurious far base; 111 = contiguous anchor base, 112..119 = the 8bp near exon.
        TestGenome genome = new TestGenome().with(CHR1, 200, 'T').set(CHR1, 101, "C").set(CHR1, 111, "CACGTACGT");

        TerminalCollapseResult res = collapser(genome).tryCollapse(CHR1, 101, "1M10N8M", bases("CACGTACGT"));

        assertTrue(res.collapsed());
        assertEquals("9M", res.newCigar());
        assertEquals(111, res.newStart());
    }

    @Test
    public void testSubAnchorJunctionClipped()
    {
        // 2M80N149M where the 2bp anchor matches neither the far exon nor the contiguous near exon (split and
        // contiguous both ~0). Contiguous ties, so the sub-threshold intron is dropped and the anchor clipped
        // -> 2S149M, start moved to the near exon. pos 81,82 intron bases do not match the read's "AC".
        TestGenome genome = new TestGenome().with(CHR1, 400, 'T').set(CHR1, 81, "GG");

        TerminalCollapseResult res = collapser(genome).tryCollapse(CHR1, 1, "2M80N149M", bases("AC" + "A".repeat(149)));

        assertTrue(res.collapsed());
        assertEquals("2S149M", res.newCigar());
        assertEquals(83, res.newStart());   // nearStart = 1 + 2 + 80
    }

    @Test
    public void testRealShortJunctionKeptWhenSplitScoresBetter()
    {
        // 50M200N5M: the 5bp anchor matches the far exon (bwa's placement, 351..355) exactly but does NOT
        // continue the near exon contiguously (151..155 mismatch). Split (+5) beats contiguous (0), so it is a
        // real short-overhang junction and is kept. near-contiguous 151..155 stay 'T'.
        TestGenome genome = new TestGenome().with(CHR1, 400, 'T').set(CHR1, 351, "ACGAC");   // far exon matches anchor

        TerminalCollapseResult res = collapser(genome)
                .tryCollapse(CHR1, 101, "50M200N5M", bases("G".repeat(50) + "ACGAC"));

        assertFalse(res.collapsed());
    }

    @Test
    public void testCleanFarAnchorWithJunkClipKept()
    {
        // exp8 shape: 50M200N5M3S. The 5M matches the far exon cleanly (split +5); the 3S matches neither the
        // far continuation nor the contiguous intron (contiguous 0). Split wins -> junction kept, clip stays.
        TestGenome genome = new TestGenome().with(CHR1, 400, 'T').set(CHR1, 351, "ACGAC");   // far exon matches the 5M anchor

        TerminalCollapseResult res = collapser(genome)
                .tryCollapse(CHR1, 101, "50M200N5M3S", bases("G".repeat(50) + "ACGAC" + "GGG"));   // 5M anchor; last 3 = junk clip

        assertFalse(res.collapsed());
    }

    @Test
    public void testAnchorAtThresholdKept()
    {
        // both flanks meet MIN_JUNCTION_ANCHOR (8): a trusted junction, neither end is a collapse candidate.
        TestGenome genome = new TestGenome().with(CHR1, 200, 'C');

        TerminalCollapseResult res = collapser(genome).tryCollapse(CHR1, 101, "8M10N8M", bases("C".repeat(16)));

        assertFalse(res.collapsed());
    }

    @Test
    public void testTrailingNonContiguousClipped()
    {
        // Anchor base doesn't continue the near exon (pos 106 is 'A', read has 'C'): reclaim nothing, drop
        // the intron, clip the anchor -> 5M1S. far base at 117 matches; contiguous pos 106 stays 'A'.
        TestGenome genome = new TestGenome().with(CHR1, 200, 'A').set(CHR1, 117, "C");

        TerminalCollapseResult res = collapser(genome).tryCollapse(CHR1, 101, "5M10N1M", bases("AAAAAC"));

        assertTrue(res.collapsed());
        assertEquals("5M1S", res.newCigar());
        assertEquals(101, res.newStart());
    }

    @Test
    public void testTrailingSoftclipFullyReclaimed()
    {
        // 5M10N1M4S: anchor + all 4 softclip bases match contiguous genome → 10M. genome 101..110.
        TestGenome genome = new TestGenome().with(CHR1, 200, 'T').set(CHR1, 101, "ACGTACGTAC");

        TerminalCollapseResult res = collapser(genome).tryCollapse(CHR1, 101, "5M10N1M4S", bases("ACGTACGTAC"));

        assertTrue(res.collapsed());
        assertEquals("10M", res.newCigar());
        assertEquals(101, res.newStart());
    }

    @Test
    public void testTrailingSoftclipPartiallyReclaimed()
    {
        // Anchor + first 3 softclip bases continue the near exon, last 2 diverge with no recovery; score
        // peaks at 3 reclaimed past the 5M near exon -> 8M2S. pos 109,110 diverge.
        TestGenome genome = new TestGenome().with(CHR1, 200, 'T').set(CHR1, 101, "ACGTACGTGG");

        TerminalCollapseResult res = collapser(genome).tryCollapse(CHR1, 101, "5M10N1M4S", bases("ACGTACGTAC"));

        assertTrue(res.collapsed());
        assertEquals("8M2S", res.newCigar());
        assertEquals(101, res.newStart());
    }

    @Test
    public void testTrailingSoftclipAnchorMismatchClipped()
    {
        // Anchor base mismatches the near-exon continuation (pos 106 = 'T'): score never goes positive, so
        // the whole terminal region is clipped and the intron dropped -> 5M5S. 101..105 set, 106 stays 'T'.
        TestGenome genome = new TestGenome().with(CHR1, 200, 'T').set(CHR1, 101, "ACGTA");

        TerminalCollapseResult res = collapser(genome).tryCollapse(CHR1, 101, "5M10N1M4S", bases("ACGTACGTAC"));

        assertTrue(res.collapsed());
        assertEquals("5M5S", res.newCigar());
        assertEquals(101, res.newStart());
    }

    @Test
    public void testLeadingSoftclipReclaimed()
    {
        // 4S1M10N8M: anchor + 4 softclip bases continue the 8bp near exon backwards contiguously; the near
        // exon (>= threshold) isn't itself a collapse candidate, so only the leading tail folds -> 13M, start 107.
        // 107..111: softclip(4) + anchor(1); 112..119: the 8bp near exon.
        TestGenome genome = new TestGenome().with(CHR1, 200, 'T').set(CHR1, 107, "CGTAC").set(CHR1, 112, "ACGTACGT");

        TerminalCollapseResult res = collapser(genome)
                .tryCollapse(CHR1, 101, "4S1M10N8M", bases("CGTACACGTACGT"));   // anchor + 8M near exon

        assertTrue(res.collapsed());
        assertEquals("13M", res.newCigar());
        assertEquals(107, res.newStart());
    }
}
